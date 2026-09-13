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
$work = Join-Path $env:TEMP ("fuelsim_h20_38_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
$job = "h20_38_hex20_sts_friction_objectivity"
$inputFile = $job + ".inp"
$resultFile = $job + ".csv"
Copy-Item (Join-Path $SourceDirectory $inputFile) $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_38.py") $work
Set-Location $work

Invoke-Abaqus `
    "job=$job" `
    "input=$inputFile" `
    output_precision=full `
    cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus H20.38 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus H20.38 solve did not report successful completion"
}

Invoke-Abaqus python extract_h20_38.py ($job + ".odb") $resultFile
if ($LASTEXITCODE -ne 0 -or !(Test-Path $resultFile)) {
    throw "Abaqus H20.38 extraction failed"
}
Copy-Item $resultFile $SourceDirectory
Write-Output "Abaqus H20.38 work directory: $work"
