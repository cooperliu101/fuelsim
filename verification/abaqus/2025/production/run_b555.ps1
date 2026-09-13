param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$JobName = "b555_small_c3d20t_friction"
$MeshName = "b549_small_c3d20t_mesh"
$Work = Join-Path $env:TEMP ("fuelsim_b555_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @(
    "$JobName.inp",
    "${MeshName}.inc",
    "${MeshName}.json",
    "extract_b555.py"
)) {
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
}
Set-Location $Work

$Stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
Invoke-Abaqus `
    job=$JobName `
    input="$JobName.inp" `
    `
    output_precision=full `
    ask_delete=OFF `
    cpus=1 interactive
$SolveExitCode = $LASTEXITCODE
$Stopwatch.Stop()
if ($SolveExitCode -ne 0) {
    throw "Abaqus B5.55 solve failed with exit code $SolveExitCode"
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
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.55 did not report successful completion"
}

Invoke-Abaqus python extract_b555.py `
    "$JobName.odb" `
    "${MeshName}.json" `
    "${JobName}_temperature.csv" `
    "${JobName}_displacement.csv" `
    "${JobName}_material.csv" `
    "${JobName}.csv"
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.55 extraction failed with exit code $LASTEXITCODE"
}

foreach ($Output in @(
    "${JobName}_timing.txt",
    "${JobName}_temperature.csv",
    "${JobName}_displacement.csv",
    "${JobName}_material.csv",
    "${JobName}.csv"
)) {
    if (!(Test-Path $Output)) {
        throw "Abaqus B5.55 extraction did not create $Output"
    }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B5.55 wall time: $($Stopwatch.Elapsed.TotalSeconds) seconds"
Write-Output "Abaqus B5.55 work directory: $Work"
