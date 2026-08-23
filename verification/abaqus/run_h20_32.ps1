param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_32_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "h20_23_hex20_surface_friction.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_32.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=h20_32_hex20_surface_friction `
    input=h20_23_hex20_surface_friction.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.32 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_32_hex20_surface_friction.sta") -or
    !(Select-String -Path "h20_32_hex20_surface_friction.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.32 solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_32.py `
    h20_32_hex20_surface_friction.odb `
    h20_32_hex20_surface_friction_contact.csv
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.32 extraction failed with exit code $LASTEXITCODE"
}

Copy-Item h20_32_hex20_surface_friction_contact.csv $SourceDirectory
Write-Output "Abaqus H20.32 work directory: $work"
