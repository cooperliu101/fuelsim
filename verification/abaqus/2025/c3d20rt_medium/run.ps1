param(
    [string]$SourceDirectory = $PSScriptRoot,
    [string]$ExistingWork = "",
    [switch]$Timing
)

$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$Job = "c3d20rt_medium_friction"
$LogPrefix = "abaqus"
if ($Timing) {
    $Job += "_timing"
    $LogPrefix += "_timing"
    # Children inherit the one-processor affinity of this PowerShell process.
    (Get-Process -Id $PID).ProcessorAffinity = [IntPtr]1
}
$Work = $ExistingWork
if (!$Work) {
    $Work = Join-Path $env:TEMP ("fuelsim_c3d20rt_medium_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
}
foreach ($Name in @("$Job.inp", "c3d20rt_medium_mesh.inc")) {
    if (!$ExistingWork) { Copy-Item (Join-Path $SourceDirectory $Name) $Work }
    if ((Get-FileHash (Join-Path $SourceDirectory $Name) -Algorithm SHA256).Hash -ne
        (Get-FileHash (Join-Path $Work $Name) -Algorithm SHA256).Hash) {
        throw "Staged file differs: $Name"
    }
}
Copy-Item (Join-Path $SourceDirectory "extract.py") $Work
Copy-Item (Join-Path $SourceDirectory "../../b548_m58_c3d20t_integrated_mesh.json") (Join-Path $Work "mesh.json")
Write-Output "work_directory=$Work"
Push-Location $Work
try {
    if (!$ExistingWork) {
        $Timer = [Diagnostics.Stopwatch]::StartNew()
        & "C:\SIMULIA\Commands\abq2025.bat" job=$Job input="$Job.inp" cpus=1 output_precision=full ask_delete=OFF interactive |
            Tee-Object -FilePath "${LogPrefix}_run.txt"
        $Code = $LASTEXITCODE
        $Timer.Stop()
        @("external_wall_seconds=$($Timer.Elapsed.TotalSeconds)", "exit_status=$Code", "cpus=1",
            "work_directory=$Work", "completed_utc=$([DateTime]::UtcNow.ToString('o'))") |
            Set-Content (Join-Path $SourceDirectory "${LogPrefix}_time.txt")
        foreach ($Name in @("$Job.dat", "$Job.msg", "$Job.sta", "${LogPrefix}_run.txt")) {
            if (Test-Path $Name) { Copy-Item $Name $SourceDirectory }
        }
        if ($Code -ne 0) { throw "Abaqus analysis returned $Code" }
    }
    if (!(Select-String -Path "$Job.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus analysis did not complete successfully"
    }
    if ($Timing) { return }
    & "C:\SIMULIA\Commands\abq2025.bat" python extract.py "$Job.odb" mesh.json reference |
        Tee-Object -FilePath "abaqus_extract.txt"
    if ($LASTEXITCODE -ne 0) { throw "Abaqus extraction failed" }
    $Frames = @(Import-Csv "reference_frames.csv")
    if ($Frames.Count -ne 20 -or [double]$Frames[-1].time -ne 1) {
        throw "Abaqus extraction did not produce all twenty frames"
    }
    foreach ($Name in @("reference_nodes.csv.gz", "reference_contact.csv.gz", "reference_points.csv.gz",
                        "reference_frames.csv", "abaqus_extract.txt")) {
        Copy-Item $Name $SourceDirectory
    }
} finally {
    Pop-Location
}
