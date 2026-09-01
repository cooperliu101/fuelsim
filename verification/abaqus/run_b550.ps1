param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b550_small_c3d20t_contact"
$MeshName = "b549_small_c3d20t_mesh"
$Work = Join-Path $env:TEMP ("fuelsim_b550_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @(
    "$JobName.inp",
    "${MeshName}.inc",
    "${MeshName}.json",
    "extract_b550.py"
)) {
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
}
Set-Location $Work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=$JobName `
    input="$JobName.inp" `
    cpus=1 `
    output_precision=full `
    ask_delete=OFF `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.50 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.50 did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" python extract_b550.py `
    "$JobName.odb" `
    "${MeshName}.json" `
    "${JobName}_temperature.csv" `
    "${JobName}_displacement.csv" `
    "${JobName}_material.csv" `
    "${JobName}.csv"
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.50 extraction failed with exit code $LASTEXITCODE"
}

foreach ($Output in @(
    "${JobName}_temperature.csv",
    "${JobName}_displacement.csv",
    "${JobName}_material.csv",
    "${JobName}.csv"
)) {
    if (!(Test-Path $Output)) {
        throw "Abaqus B5.50 extraction did not create $Output"
    }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B5.50 work directory: $Work"
