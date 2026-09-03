param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b551_m58_c3d20t_finite_sliding"
$MeshName = "b548_m58_c3d20t_integrated"
$Work = Join-Path $env:TEMP ("fuelsim_b551_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @(
    "$JobName.inp",
    "${MeshName}_mesh.inc",
    "${MeshName}_mesh.json",
    "extract_b548.py"
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
    throw "Abaqus B5.51 solve failed with exit code $SolveExitCode"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.51 did not report successful completion"
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

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b548.py `
    "$JobName.odb" `
    "${MeshName}_mesh.json" `
    "${JobName}_nodal.csv" `
    "${JobName}_contact.csv" `
    "${JobName}_clad_points.csv"
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.51 field extraction failed with exit code $LASTEXITCODE"
}

foreach ($Output in @(
    "${JobName}_nodal.csv",
    "${JobName}_contact.csv",
    "${JobName}_clad_points.csv"
)) {
    if (!(Test-Path $Output)) {
        throw "Abaqus B5.51 field extraction did not create $Output"
    }
}

foreach ($Output in @(
    "${JobName}_timing.txt",
    "${JobName}_nodal.csv",
    "${JobName}_contact.csv",
    "${JobName}_clad_points.csv"
)) {
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B5.51 wall time: $($Stopwatch.Elapsed.TotalSeconds) seconds"
Write-Output "Abaqus B5.51 work directory: $Work"
