param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b548_m58_c3d20t_integrated"
$Work = Join-Path $env:TEMP ("fuelsim_b548_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @(
    "$JobName.inp",
    "${JobName}_mesh.inc"
)) {
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
}
Set-Location $Work

$Stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
& "C:\SIMULIA\Commands\abaqus.bat" `
    job=$JobName `
    input="$JobName.inp" `
    cpus=1 `
    output_precision=full `
    ask_delete=OFF `
    interactive
$SolveExitCode = $LASTEXITCODE
$Stopwatch.Stop()
if ($SolveExitCode -ne 0) {
    throw "Abaqus B5.48 solve failed with exit code $SolveExitCode"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.48 did not report successful completion"
}

$Timing = @(
    "abaqus_external_wall_seconds=$($Stopwatch.Elapsed.TotalSeconds.ToString('F6', [System.Globalization.CultureInfo]::InvariantCulture))",
    "abaqus_cpus=1",
    "abaqus_output_precision=full",
    "windows_computer=$env:COMPUTERNAME"
)
if (Test-Path "$JobName.dat") {
    $DatLines = Get-Content "$JobName.dat"
    $SummaryStart = -1
    for ($Index = 0; $Index -lt $DatLines.Count; ++$Index) {
        if ($DatLines[$Index] -match "JOB TIME SUMMARY") {
            $SummaryStart = $Index
        }
    }
    if ($SummaryStart -ge 0) {
        $Timing += "abaqus_job_time_summary_begin"
        $SummaryEnd = [Math]::Min($SummaryStart + 4, $DatLines.Count - 1)
        $Timing += $DatLines[$SummaryStart..$SummaryEnd]
        $Timing += "abaqus_job_time_summary_end"
    }
}
$Timing | Set-Content -Encoding ASCII "${JobName}_timing.txt"

$IncrementSummary = @(
    "increment`tattempt`tsevere_discontinuity_iterations`tequilibrium_iterations`ttotal_iterations`ttime_s`tincrement_s"
)
foreach ($Line in Get-Content "$JobName.sta") {
    if ($Line -match "^\s+1\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+") {
        $Fields = ($Line.Trim() -split "\s+")
        $IncrementSummary += "$($Fields[1])`t$($Fields[2])`t$($Fields[3])`t$($Fields[4])`t$($Fields[5])`t$($Fields[7])`t$($Fields[8])"
    }
}
if ($IncrementSummary.Count -ne 21) {
    throw "Abaqus B5.48 did not produce twenty increment summary rows"
}
$IncrementSummary | Set-Content -Encoding ASCII "${JobName}_increments.tsv"

foreach ($Output in @("${JobName}_timing.txt", "${JobName}_increments.tsv")) {
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B5.48 wall time: $($Stopwatch.Elapsed.TotalSeconds) seconds"
Write-Output "Abaqus B5.48 work directory: $Work"
