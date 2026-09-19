param(
    [Parameter(Mandatory=$true)][string]$SourceDirectory,
    [Parameter(Mandatory=$true)][string]$DestinationDirectory,
    [Parameter(Mandatory=$true)][string]$CaseName,
    [switch]$Extract
)
$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:NUMEXPR_NUM_THREADS = "1"
$work = Join-Path $env:TEMP ("fuelsim_thermal_perf_"+[guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Force -Path $work,$DestinationDirectory | Out-Null
Copy-Item (Join-Path $SourceDirectory "*.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract.py") $work
$arguments = '/d /c start "" /b /wait /affinity 1 "'+$env:ComSpec+
    '" /d /c C:\SIMULIA\Commands\abaqus.bat'+
    " job=$CaseName input=$CaseName.inp cpus=1 ask_delete=OFF interactive"
$watch = [System.Diagnostics.Stopwatch]::StartNew()
$process = Start-Process -FilePath $env:ComSpec -ArgumentList $arguments -WorkingDirectory $work -Wait -PassThru
$watch.Stop()
if ($process.ExitCode -ne 0) { throw "Abaqus exit code $($process.ExitCode)" }
if (!(Select-String -Path (Join-Path $work "$CaseName.sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus did not complete"
}
@{external_wall_seconds=$watch.Elapsed.TotalSeconds; cpus=1; affinity_mask=1;
  computer=$env:COMPUTERNAME; work_directory=$work;
  input_sha256=(Get-FileHash (Join-Path $work "$CaseName.inp") -Algorithm SHA256).Hash;
  processors=@(Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors)} |
  ConvertTo-Json -Depth 5 | Set-Content -Encoding ASCII (Join-Path $DestinationDirectory "timing.json")
foreach ($extension in @("dat","sta","msg","log")) {
    if (Test-Path (Join-Path $work "$CaseName.$extension")) {
        Copy-Item (Join-Path $work "$CaseName.$extension") $DestinationDirectory
    }
}
if ($Extract) {
    Push-Location $work
    try {
        & "C:\SIMULIA\Commands\abaqus.bat" python extract.py "$CaseName.odb" "$CaseName.npz"
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed" }
        Copy-Item "$CaseName.npz" $DestinationDirectory
    } finally { Pop-Location }
}
Write-Output "case=$CaseName external_wall_seconds=$($watch.Elapsed.TotalSeconds) work=$work"
