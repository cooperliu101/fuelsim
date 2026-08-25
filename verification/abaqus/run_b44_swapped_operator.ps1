param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b44_operator_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
$job = "b44_hex8_swapped_operator_probe"
$inputFile = $job + ".inp"
$outputFile = $job + ".csv"
Copy-Item (Join-Path $SourceDirectory $inputFile) $work
Copy-Item (Join-Path $SourceDirectory "extract_b44_swapped_operator.py") $work
Set-Location $work
& "C:\SIMULIA\Commands\abaqus.bat" `
    "job=$job" `
    "input=$inputFile" `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.4 operator solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.4 operator probe did not report successful completion"
}
& "C:\SIMULIA\Commands\abaqus.bat" python extract_b44_swapped_operator.py ($job + ".odb") $outputFile
if ($LASTEXITCODE -ne 0 -or !(Test-Path $outputFile) -or (Get-Content $outputFile).Count -ne 280) {
    throw "Abaqus B4.4 operator extraction failed"
}
Copy-Item $outputFile $SourceDirectory
Write-Output "Abaqus B4.4 operator work directory: $work"
