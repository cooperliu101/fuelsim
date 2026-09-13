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
$work = Join-Path $env:TEMP ("fuelsim_b47_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_b47.py") $work
Set-Location $work

foreach ($name in @("unit", "aspect16", "area4", "area_quarter")) {
    $job = "b47_hex8_scale_" + $name
    $inputFile = $job + ".inp"
    $nodeFile = $job + "_nodes.csv"
    $contactFile = $job + "_contact.csv"
    Copy-Item (Join-Path $SourceDirectory $inputFile) $work
    Invoke-Abaqus `
        "job=$job" `
        "input=$inputFile" `
        output_precision=full `
        cpus=1 interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus B4.7 $name solve failed with exit code $LASTEXITCODE"
    }
    if (!(Test-Path ($job + ".sta")) -or
        !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus B4.7 $name did not report successful completion"
    }
    Invoke-Abaqus python extract_b47.py ($job + ".odb") $nodeFile $contactFile
    if ($LASTEXITCODE -ne 0 -or !(Test-Path $nodeFile) -or !(Test-Path $contactFile) -or
        (Get-Content $nodeFile).Count -ne 49 -or (Get-Content $contactFile).Count -ne 13) {
        throw "Abaqus B4.7 $name extraction failed"
    }
    Copy-Item $nodeFile $SourceDirectory
    Copy-Item $contactFile $SourceDirectory
}
Write-Output "Abaqus B4.7 work directory: $work"
