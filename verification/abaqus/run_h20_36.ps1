param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_36_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
$job = "h20_36_hex20_sts_quadratic_cylinder_friction_path"
$inputFile = $job + ".inp"
Copy-Item (Join-Path $SourceDirectory $inputFile) $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_36.py") $work
Set-Location $work

& "C:\SIMULIA\Commands\abaqus.bat" `
    "job=$job" `
    "input=$inputFile" `
    output_precision=full `
    interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.36 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.36 solve did not report successful completion"
}

& "C:\SIMULIA\Commands\abaqus.bat" `
    python extract_h20_36.py `
    ($job + ".odb") `
    $inputFile `
    ($job + "_displacement.csv") `
    ($job + "_contact.csv")
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.36 extraction failed with exit code $LASTEXITCODE"
}

foreach ($result in @(($job + "_displacement.csv"), ($job + "_contact.csv"))) {
    if (!(Test-Path $result)) {
        throw "Abaqus H20.36 extraction did not create $result"
    }
    Copy-Item $result $SourceDirectory
}
Write-Output "Abaqus H20.36 work directory: $work"
