param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b43_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
$job = "b43_hex8_finite_strain_contact"
$inputFile = $job + ".inp"
$nodeFile = $job + "_nodes.csv"
$contactFile = $job + ".csv"
Copy-Item (Join-Path $SourceDirectory $inputFile) $work
Copy-Item (Join-Path $SourceDirectory "extract_b43.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    "job=$job" `
    "input=$inputFile" `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.3 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.3 solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b43.py ($job + ".odb") $nodeFile $contactFile
if ($LASTEXITCODE -ne 0 -or !(Test-Path $nodeFile) -or !(Test-Path $contactFile) -or
    (Get-Content $nodeFile).Count -ne 321 -or (Get-Content $contactFile).Count -ne 81) {
    throw "Abaqus B4.3 extraction failed"
}
Copy-Item $nodeFile $SourceDirectory
Copy-Item $contactFile $SourceDirectory
Write-Output "Abaqus B4.3 work directory: $work"
