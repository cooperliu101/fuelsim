param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b46_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_b46.py") $work
$job = "b46_hex8_release_recontact"
$inputFile = $job + ".inp"
$nodeFile = $job + "_nodes.csv"
$contactFile = $job + "_contact.csv"
Copy-Item (Join-Path $SourceDirectory $inputFile) $work
Set-Location $work
& "C:\SIMULIA\Commands\abaqus.bat" `
    "job=$job" `
    "input=$inputFile" `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.6 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.6 did not report successful completion"
}
& "C:\SIMULIA\Commands\abaqus.bat" python extract_b46.py ($job + ".odb") $nodeFile $contactFile
if ($LASTEXITCODE -ne 0 -or !(Test-Path $nodeFile) -or !(Test-Path $contactFile) -or
    (Get-Content $nodeFile).Count -ne 121 -or (Get-Content $contactFile).Count -ne 25) {
    throw "Abaqus B4.6 extraction failed"
}
Copy-Item $nodeFile $SourceDirectory
Copy-Item $contactFile $SourceDirectory
Write-Output "Abaqus B4.6 work directory: $work"
