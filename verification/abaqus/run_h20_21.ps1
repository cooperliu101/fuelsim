param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_21_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "h20_21_hex20_surface_contact.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_21.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=h20_21_hex20_surface_contact `
    input=h20_21_hex20_surface_contact.inp `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_21_hex20_surface_contact.sta") -or
    !(Select-String -Path "h20_21_hex20_surface_contact.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_21.py `
    h20_21_hex20_surface_contact.odb `
    h20_21_hex20_surface_contact_displacement.csv `
    h20_21_hex20_surface_contact_reaction.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus extraction failed with exit code $LASTEXITCODE"
}

Copy-Item h20_21_hex20_surface_contact_displacement.csv $SourceDirectory
Copy-Item h20_21_hex20_surface_contact_reaction.csv $SourceDirectory
Write-Output "Abaqus H20.21 work directory: $work"
