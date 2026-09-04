param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory,
    [ValidateSet("steady", "transient", "ramped")]
    [string]$Mode = "ramped"
)

$ErrorActionPreference = "Stop"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
switch ($Mode) {
    "steady" {
        $Stem = "b60_fuel_plate_c3d20t_finite_steady_bending"
        $FinalTime = "1.0"
    }
    "transient" {
        $Stem = "b60_fuel_plate_c3d20t_finite_bending"
        $FinalTime = "10.0"
    }
    "ramped" {
        $Stem = "b60_fuel_plate_c3d20t_finite_ramped_bending"
        $FinalTime = "10.0"
    }
}
$Work = Join-Path $env:TEMP ("fuelsim_" + $Stem + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @(
    "$Stem.inp",
    "b60_long_plate_meat_clad_c3d20t_mesh.inc",
    "b60_long_plate_meat_clad_c3d20t_mesh.json",
    "extract_b60_c3d20t.py"
)) {
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
}
Set-Location $Work

$Stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
& "C:\SIMULIA\Commands\abaqus.bat" `
    job=$Stem `
    input="$Stem.inp" `
    cpus=1 `
    output_precision=full `
    ask_delete=OFF `
    interactive
$SolveExitCode = $LASTEXITCODE
$Stopwatch.Stop()
if ($SolveExitCode -ne 0) { throw "Abaqus B6.0 C3D20T $Mode solve failed with exit code $SolveExitCode" }
if (!(Test-Path "$Stem.sta") -or
    !(Select-String -Path "$Stem.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B6.0 C3D20T $Mode solve did not report successful completion"
}

$Timing = @(
    "abaqus_external_wall_seconds=$($Stopwatch.Elapsed.TotalSeconds.ToString('F6', [System.Globalization.CultureInfo]::InvariantCulture))",
    "abaqus_cpus=1",
    "abaqus_output_precision=full",
    "windows_computer=$env:COMPUTERNAME"
)
if (Test-Path "$Stem.dat") {
    $DatLines = Get-Content "$Stem.dat"
    $SummaryStart = -1
    for ($Index = 0; $Index -lt $DatLines.Count; ++$Index) {
        if ($DatLines[$Index] -match "JOB TIME SUMMARY") { $SummaryStart = $Index }
    }
    if ($SummaryStart -ge 0) {
        $Timing += "abaqus_job_time_summary_begin"
        $SummaryEnd = [Math]::Min($SummaryStart + 4, $DatLines.Count - 1)
        $Timing += @($DatLines[$SummaryStart..$SummaryEnd] | ForEach-Object { $_.TrimEnd() })
        $Timing += "abaqus_job_time_summary_end"
    }
}
[System.IO.File]::WriteAllText("$PWD\${Stem}_timing.txt", (($Timing -join "`n") + "`n"), [System.Text.Encoding]::ASCII)

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b60_c3d20t.py `
    "$Stem.odb" "b60_long_plate_meat_clad_c3d20t_mesh.json" `
    "${Stem}_temperature.csv" "${Stem}_displacement.csv" "${Stem}_material.csv"
if ($LASTEXITCODE -ne 0) { throw "Abaqus B6.0 C3D20T $Mode extraction failed with exit code $LASTEXITCODE" }
foreach ($Output in @(
    "${Stem}_temperature.csv",
    "${Stem}_displacement.csv",
    "${Stem}_material.csv",
    "${Stem}_timing.txt"
)) {
    if (!(Test-Path $Output)) { throw "Abaqus B6.0 C3D20T $Mode did not create $Output" }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B6.0 C3D20T $Mode wall time: $($Stopwatch.Elapsed.TotalSeconds) seconds"
Write-Output "Abaqus B6.0 C3D20T $Mode work directory: $Work"
