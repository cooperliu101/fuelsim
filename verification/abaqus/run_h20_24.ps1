param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_24_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "h20_24_hex20_nonmatching_surface_contact.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_24.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=h20_24_hex20_nonmatching_surface_contact `
    input=h20_24_hex20_nonmatching_surface_contact.inp `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_24_hex20_nonmatching_surface_contact.sta") -or
    !(Select-String -Path "h20_24_hex20_nonmatching_surface_contact.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_24.py `
    h20_24_hex20_nonmatching_surface_contact.odb `
    h20_24_hex20_nonmatching_surface_contact_displacement.csv `
    h20_24_hex20_nonmatching_surface_contact_reaction.csv `
    h20_24_hex20_nonmatching_surface_contact_pressure.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus extraction failed with exit code $LASTEXITCODE"
}
foreach ($result in @(
    "h20_24_hex20_nonmatching_surface_contact_displacement.csv",
    "h20_24_hex20_nonmatching_surface_contact_reaction.csv",
    "h20_24_hex20_nonmatching_surface_contact_pressure.csv")) {
    if (!(Test-Path $result)) {
        throw "Abaqus extraction did not create $result"
    }
}

Copy-Item h20_24_hex20_nonmatching_surface_contact_displacement.csv $SourceDirectory
Copy-Item h20_24_hex20_nonmatching_surface_contact_reaction.csv $SourceDirectory
Copy-Item h20_24_hex20_nonmatching_surface_contact_pressure.csv $SourceDirectory
Write-Output "Abaqus H20.24 work directory: $work"
