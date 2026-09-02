param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "h20_41_hex20_finite_sliding_partial_contact"
$Work = Join-Path $env:TEMP ("fuelsim_h20_41_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @("$JobName.inp", "extract_h20_41.py")) {
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
}
Set-Location $Work

& "C:\SIMULIA\Commands\abaqus.bat" `
    "job=$JobName" `
    "input=$JobName.inp" `
    cpus=1 `
    output_precision=full `
    ask_delete=OFF `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.41 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.41 did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" python extract_h20_41.py `
    "$JobName.odb" `
    "$JobName.inp" `
    "${JobName}_contact.csv"
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.41 extraction failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${JobName}_contact.csv")) {
    throw "Abaqus H20.41 extraction did not create the contact reference"
}
Copy-Item "${JobName}_contact.csv" $SourceDirectory
Write-Output "Abaqus H20.41 work directory: $Work"
