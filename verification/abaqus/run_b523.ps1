param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b523_hex8_c3d8t_integrated_path"
$Work = Join-Path $env:TEMP ("fuelsim_b523_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_b523.py") $Work
Set-Location $Work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=$JobName `
    input="$JobName.inp" `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.23 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.23 did not report successful completion"
}

$Outputs = @(
    "${JobName}_nodal.csv",
    "${JobName}_integration.csv",
    "${JobName}_contact.csv",
    "${JobName}_energy.csv"
)
& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_b523.py `
    "$JobName.odb" `
    $Outputs[0] `
    $Outputs[1] `
    $Outputs[2] `
    $Outputs[3]
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.23 extraction failed with exit code $LASTEXITCODE"
}
foreach ($Output in $Outputs) {
    if (!(Test-Path $Output)) {
        throw "Abaqus B5.23 extraction did not create $Output"
    }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B5.23 work directory: $Work"
