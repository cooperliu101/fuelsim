param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_33_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "h20_33_hex20_nonmatching_friction_path.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_33.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=h20_33_hex20_nonmatching_friction_path `
    input=h20_33_hex20_nonmatching_friction_path.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.33 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_33_hex20_nonmatching_friction_path.sta") -or
    !(Select-String -Path "h20_33_hex20_nonmatching_friction_path.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.33 solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_33.py `
    h20_33_hex20_nonmatching_friction_path.odb `
    h20_33_hex20_nonmatching_friction_path_displacement.csv `
    h20_33_hex20_nonmatching_friction_path_contact.csv `
    h20_33_hex20_nonmatching_friction_path_reaction.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.33 extraction failed with exit code $LASTEXITCODE"
}

foreach ($result in @(
    "h20_33_hex20_nonmatching_friction_path_displacement.csv",
    "h20_33_hex20_nonmatching_friction_path_contact.csv",
    "h20_33_hex20_nonmatching_friction_path_reaction.csv")) {
    if (!(Test-Path $result)) {
        throw "Abaqus H20.33 extraction did not create $result"
    }
    Copy-Item $result $SourceDirectory
}
Write-Output "Abaqus H20.33 work directory: $work"
