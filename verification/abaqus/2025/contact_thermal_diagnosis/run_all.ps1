param([string]$TaskFile = (Join-Path $PSScriptRoot "tasks.json"), [int]$MaximumJobs = 8, [string]$LogDirectory = "refresh_logs")
$ErrorActionPreference = "Stop"
if ($MaximumJobs -lt 1 -or $MaximumJobs -gt 8) { throw "Native task concurrency must be between one and eight" }
$Root = $PSScriptRoot
$Logs = Join-Path $Root $LogDirectory
New-Item -ItemType Directory -Force -Path $Logs | Out-Null
$Tasks = Get-Content $TaskFile -Raw | ConvertFrom-Json
$Active = @(); $Results = @(); $Index = 0
while ($Index -lt $Tasks.Count -or $Active.Count -gt 0) {
    while ($Index -lt $Tasks.Count -and $Active.Count -lt $MaximumJobs) {
        $Task = $Tasks[$Index]; $Index++
        $Active += Start-Job -Name $Task.id -ArgumentList $Root, $Logs, $Task -ScriptBlock {
            param($Root, $Logs, $Task)
            $ErrorActionPreference = "Stop"
            $env:OMP_NUM_THREADS = "1"; $env:MKL_NUM_THREADS = "1"; $env:OPENBLAS_NUM_THREADS = "1"
            $Source = if ($Task.directory) { Join-Path $Root $Task.directory } else { $Root }
            $Arguments = @{}
            foreach ($Property in $Task.args.PSObject.Properties) {
                $Arguments[$Property.Name] = if ($Property.Value -eq '$SOURCE') { $Source } else { $Property.Value }
            }
            $Log = Join-Path $Logs ($Task.id + ".log")
            $Started = [DateTime]::UtcNow.ToString('o')
            try {
                & (Join-Path $Source $Task.script) @Arguments *>&1 | Out-File -FilePath $Log -Encoding utf8
                $Body = if (Test-Path $Log) { Get-Content $Log -Raw } else { "" }
                if ($Body -match 'Traceback \(most recent call last\)|Abaqus Error:|ANALYSIS HAS NOT BEEN COMPLETED') { throw "Native runner log contains a solve or extraction failure" }
                [pscustomobject]@{task=$Task.id; status="success"; started=$Started; finished=[DateTime]::UtcNow.ToString('o'); error=""}
            } catch {
                $_ | Out-File -FilePath $Log -Append -Encoding utf8
                [pscustomobject]@{task=$Task.id; status="failed"; started=$Started; finished=[DateTime]::UtcNow.ToString('o'); error=$_.ToString()}
            }
        }
    }
    $Finished = @($Active | Where-Object { $_.State -in @('Completed', 'Failed', 'Stopped') })
    foreach ($Job in $Finished) {
        $Record = Receive-Job $Job
        if (!$Record) { $Record = [pscustomobject]@{task=$Job.Name; status="failed"; error="Worker returned no status"} }
        $Results += $Record
        Write-Output ("{0}/{1} {2} {3}" -f $Results.Count, $Tasks.Count, $Job.Name, $Record.status)
        Remove-Job $Job
    }
    $FinishedIds = @($Finished | ForEach-Object { $_.Id })
    $Active = @($Active | Where-Object { $_.Id -notin $FinishedIds })
    $Results | Export-Csv (Join-Path $Logs "status.csv") -NoTypeInformation -Encoding utf8
    if ($Active.Count -gt 0) { Wait-Job -Job $Active -Any -Timeout 1 | Out-Null }
}

if ($Results.Count -ne $Tasks.Count -or @($Results | Where-Object { $_.status -ne "success" }).Count -gt 0) { throw "One or more native tasks failed; inspect the status CSV and logs" }
