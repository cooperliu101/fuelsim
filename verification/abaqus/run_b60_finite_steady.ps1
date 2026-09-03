param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b60_fuel_plate_c3d8rt_finite_steady_bending"
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$Work = Join-Path $env:TEMP ("fuelsim_b60_steady_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @(
    "b60_fuel_plate_c3d8rt_finite_steady_bending.inp",
    "b60_long_plate_meat_clad_c3d8rt.inc",
    "b60_long_plate_meat_clad_c3d8rt.json",
    "extract_b60_finite.py"
)) {
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
}
Set-Location $Work

$Stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
& "C:\SIMULIA\Commands\abaqus.bat" `
    job=$JobName `
    input="b60_fuel_plate_c3d8rt_finite_steady_bending.inp" `
    cpus=1 `
    output_precision=full `
    ask_delete=OFF `
    interactive
$SolveExitCode = $LASTEXITCODE
$Stopwatch.Stop()
if ($SolveExitCode -ne 0) {
    throw "Abaqus B6.0 steady finite-strain solve failed with exit code $SolveExitCode"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B6.0 steady finite-strain solve did not report successful completion"
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
        if ($DatLines[$Index] -match "JOB TIME SUMMARY") { $SummaryStart = $Index }
    }
    if ($SummaryStart -ge 0) {
        $Timing += "abaqus_job_time_summary_begin"
        $SummaryEnd = [Math]::Min($SummaryStart + 4, $DatLines.Count - 1)
        $Timing += @($DatLines[$SummaryStart..$SummaryEnd] | ForEach-Object { $_.TrimEnd() })
        $Timing += "abaqus_job_time_summary_end"
    }
}
if (Test-Path "$JobName.msg") {
    $MessageText = Get-Content "$JobName.msg" -Raw
    foreach ($Metric in @(
        @{ Name = "abaqus_increments"; Pattern = "TOTAL OF\s+(\d+)\s+INCREMENTS" },
        @{ Name = "abaqus_nonlinear_iterations"; Pattern = "(\d+)\s+ITERATIONS INCLUDING CONTACT ITERATIONS" },
        @{ Name = "abaqus_equation_solver_passes"; Pattern = "(\d+)\s+PASSES THROUGH THE EQUATION SOLVER" },
        @{ Name = "abaqus_matrix_decompositions"; Pattern = "(\d+)\s+INVOLVE MATRIX DECOMPOSITION" },
        @{ Name = "abaqus_equation_reorderings"; Pattern = "(\d+)\s+REORDERING OF EQUATIONS" },
        @{ Name = "abaqus_line_search_residual_evaluations"; Pattern = "(\d+)\s+ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES" },
        @{ Name = "abaqus_line_search_operator_evaluations"; Pattern = "(\d+)\s+ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES" }
    )) {
        if ($MessageText -match $Metric.Pattern) { $Timing += "$($Metric.Name)=$($Matches[1])" }
    }
}
$TimingText = ($Timing -join "`n") + "`n"
[System.IO.File]::WriteAllText("$PWD\${JobName}_timing.txt", $TimingText, [System.Text.Encoding]::ASCII)

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b60_finite.py `
    "$JobName.odb" "b60_long_plate_meat_clad_c3d8rt.json" "${JobName}_nodal.csv" 1.0
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B6.0 steady finite-strain extraction failed with exit code $LASTEXITCODE"
}
foreach ($Output in @("${JobName}_nodal.csv", "${JobName}_timing.txt")) {
    if (!(Test-Path $Output)) { throw "Abaqus B6.0 steady finite-strain solve did not create $Output" }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B6.0 steady finite-strain wall time: $($Stopwatch.Elapsed.TotalSeconds) seconds"
Write-Output "Abaqus B6.0 steady finite-strain work directory: $Work"
