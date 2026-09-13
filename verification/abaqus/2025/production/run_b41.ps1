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
$work = Join-Path $env:TEMP ("fuelsim_b41_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_b41.py") $work
Set-Location $work

foreach ($mode in @("small", "finite")) {
    $job = "b41_hex8_finite_sliding_" + $mode
    $inputFile = $job + ".inp"
    $resultFile = $job + ".csv"
    Copy-Item (Join-Path $SourceDirectory $inputFile) $work
    Invoke-Abaqus `
        "job=$job" `
        "input=$inputFile" `
        output_precision=full `
        cpus=1 interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus B4.1 $mode solve failed with exit code $LASTEXITCODE"
    }
    if (!(Test-Path ($job + ".sta")) -or
        !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus B4.1 $mode solve did not report successful completion"
    }
    Invoke-Abaqus python extract_b41.py ($job + ".odb") $resultFile
    if ($LASTEXITCODE -ne 0 -or !(Test-Path $resultFile) -or (Get-Item $resultFile).Length -le 128) {
        throw "Abaqus B4.1 $mode extraction failed"
    }
    Copy-Item $resultFile $SourceDirectory
}
Write-Output "Abaqus B4.1 work directory: $work"
