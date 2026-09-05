param(
    [string]$SourceDirectory,
    [Parameter(Mandatory = $true)][string]$DestinationDirectory,
    [int]$Repetitions = 3
)

$ErrorActionPreference = "Stop"
if ([string]::IsNullOrWhiteSpace($SourceDirectory)) { $SourceDirectory = $PSScriptRoot }
if ($Repetitions -lt 1) { throw "Repetitions must be positive" }
$JobName = "b556_m58_c3d20t_finite_sliding_friction_timing"
$MeshName = "b548_m58_c3d20t_integrated_mesh.inc"
$env:OMP_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:NUMEXPR_NUM_THREADS = "1"
New-Item -ItemType Directory -Force -Path $DestinationDirectory | Out-Null

for ($Sample = 1; $Sample -le $Repetitions; ++$Sample) {
    $Work = Join-Path $env:TEMP ("fuelsim_b556_timing_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    foreach ($Name in @("$JobName.inp", $MeshName)) {
        Copy-Item (Join-Path $SourceDirectory $Name) $Work
    }
    Push-Location $Work
    try {
        Write-Output "abaqus_sample=$Sample begin work_directory=$Work"
        # START applies affinity when creating the child, before Abaqus can
        # launch its analysis process. WAIT includes the complete batch job.
        $Arguments = '/d /c start "" /b /wait /affinity 1 "' + $env:ComSpec +
            '" /d /c C:\SIMULIA\Commands\abaqus.bat' +
            " job=$JobName input=$JobName.inp cpus=1 output_precision=full ask_delete=OFF interactive"
        $Watch = [System.Diagnostics.Stopwatch]::StartNew()
        $Process = Start-Process -FilePath $env:ComSpec -ArgumentList $Arguments -WorkingDirectory $Work -PassThru -Wait
        $Watch.Stop()
        if ($Process.ExitCode -ne 0) { throw "Abaqus failed with exit code $($Process.ExitCode)" }
        if (!(Test-Path "$JobName.sta") -or
            !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
            throw "Abaqus did not report successful completion"
        }
        $Seconds = $Watch.Elapsed.TotalSeconds.ToString("F6", [System.Globalization.CultureInfo]::InvariantCulture)
        $Timing = @("sample=$Sample", "abaqus_external_wall_seconds=$Seconds", "cpus=1",
            "processor_affinity=1", "numerical_library_threads=1", "field_output=disabled",
            "history_output=disabled", "restart_output=disabled", "windows_computer=$env:COMPUTERNAME",
            "input_sha256=$((Get-FileHash "$JobName.inp" -Algorithm SHA256).Hash)")
        $Lines = Get-Content "$JobName.dat"
        $SummaryStart = -1
        for ($Index = 0; $Index -lt $Lines.Count; ++$Index) {
            if ($Lines[$Index] -match "JOB TIME SUMMARY") { $SummaryStart = $Index }
        }
        if ($SummaryStart -ge 0) {
            $Timing += $Lines[$SummaryStart..([Math]::Min($SummaryStart + 4, $Lines.Count - 1))]
        }
        foreach ($Extension in @("sta", "dat", "msg")) {
            Copy-Item "$JobName.$Extension" (Join-Path $DestinationDirectory "abaqus_sample_$Sample.$Extension")
        }
        if (Test-Path "$JobName.log") {
            Copy-Item "$JobName.log" (Join-Path $DestinationDirectory "abaqus_sample_$Sample.log")
        }
        $Timing | Set-Content -Encoding ASCII (Join-Path $DestinationDirectory "abaqus_sample_$Sample.txt")
        Write-Output "abaqus_sample=$Sample external_seconds=$Seconds work_directory=$Work"
    } finally {
        Pop-Location
    }
}
