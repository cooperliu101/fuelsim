param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$job = "h20_31_hex20_curved_sts_operator"
$work = Join-Path $env:TEMP ("fuelsim_h20_31_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory ($job + ".inp")) $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_31.py") $work
Set-Location $work
& "C:\SIMULIA\Commands\abaqus.bat" `
    "job=$job" `
    ("input=" + $job + ".inp") `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.31 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.31 did not report successful completion"
}
& "C:\SIMULIA\Commands\abaqus.bat" python extract_h20_31.py ($job + ".odb") ($job + ".csv")
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.31 extraction failed with exit code $LASTEXITCODE"
}
Copy-Item ($job + ".csv") $SourceDirectory
Write-Output "Abaqus H20.31 work directory: $work"
