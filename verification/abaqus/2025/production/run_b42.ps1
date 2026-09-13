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
$work = Join-Path $env:TEMP ("fuelsim_b42_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
$job = "b42_hex8_sts_friction_objectivity"
$inputFile = $job + ".inp"
$resultFile = $job + ".csv"
Copy-Item (Join-Path $SourceDirectory $inputFile) $work
Copy-Item (Join-Path $SourceDirectory "extract_b42.py") $work
Set-Location $work

Invoke-Abaqus `
    "job=$job" `
    "input=$inputFile" `
    output_precision=full `
    cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B4.2 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path ($job + ".sta")) -or
    !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B4.2 solve did not report successful completion"
}

Invoke-Abaqus python extract_b42.py ($job + ".odb") $resultFile
if ($LASTEXITCODE -ne 0 -or !(Test-Path $resultFile)) {
    throw "Abaqus B4.2 extraction failed"
}
Copy-Item $resultFile $SourceDirectory
Write-Output "Abaqus B4.2 work directory: $work"
