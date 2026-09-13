param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_35_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
$job = "h20_35_hex20_sts_quadratic_cylinder_friction"
$inputFile = $job + ".inp"
Copy-Item (Join-Path $SourceDirectory $inputFile) $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_35.py") $work
Set-Location $work

Invoke-Abaqus `
    "job=$job" `
    "input=$inputFile" `
    output_precision=full `
    cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.35 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.35 solve did not report successful completion"
}

Invoke-Abaqus `
    python extract_h20_35.py `
    ($job + ".odb") `
    $inputFile `
    ($job + "_displacement.csv") `
    ($job + "_contact.csv")
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.35 extraction failed with exit code $LASTEXITCODE"
}

foreach ($result in @(($job + "_displacement.csv"), ($job + "_contact.csv"))) {
    if (!(Test-Path $result)) {
        throw "Abaqus H20.35 extraction did not create $result"
    }
    Copy-Item $result $SourceDirectory
}
Write-Output "Abaqus H20.35 work directory: $work"
