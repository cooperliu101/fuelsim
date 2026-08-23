param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_27_transfer_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "h20_24_hex20_nonmatching_surface_contact.inp") $work
Copy-Item (Join-Path $SourceDirectory "generate_h20_27_transfer.py") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_27_transfer.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    python generate_h20_27_transfer.py `
    h20_24_hex20_nonmatching_surface_contact.inp `
    h20_27_hex20_small_sliding_transfer.inp
if (!(Test-Path "h20_27_hex20_small_sliding_transfer.inp")) {
    throw "Abaqus H20.27 transfer input generation failed"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    job=h20_27_hex20_small_sliding_transfer `
    input=h20_27_hex20_small_sliding_transfer.inp `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.27 transfer solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "h20_27_hex20_small_sliding_transfer.sta") -or
    !(Select-String -Path "h20_27_hex20_small_sliding_transfer.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.27 transfer solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_27_transfer.py `
    h20_27_hex20_small_sliding_transfer.odb `
    h20_27_hex20_small_sliding_transfer_nodal.csv `
    h20_27_hex20_small_sliding_transfer_history.csv
if (!(Test-Path "h20_27_hex20_small_sliding_transfer_nodal.csv") -or
    !(Test-Path "h20_27_hex20_small_sliding_transfer_history.csv")) {
    throw "Abaqus H20.27 transfer extraction did not create both CSV files"
}

Copy-Item h20_27_hex20_small_sliding_transfer.inp $SourceDirectory
Copy-Item h20_27_hex20_small_sliding_transfer_nodal.csv $SourceDirectory
Copy-Item h20_27_hex20_small_sliding_transfer_history.csv $SourceDirectory
Write-Output "Abaqus H20.27 transfer work directory: $work"
