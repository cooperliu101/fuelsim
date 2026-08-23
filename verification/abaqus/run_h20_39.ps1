param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_39_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_h20_39.py") $work
Set-Location $work

foreach ($mode in @("small", "finite")) {
    $job = "h20_39_hex20_finite_sliding_" + $mode
    $inputFile = $job + ".inp"
    $resultFile = $job + ".csv"
    Copy-Item (Join-Path $SourceDirectory $inputFile) $work
    & "C:\SIMULIA\Commands\abaqus.bat" `
        "job=$job" `
        "input=$inputFile" `
        output_precision=full `
        interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus H20.39 $mode solve failed with exit code $LASTEXITCODE"
    }
    if (!(Test-Path ($job + ".sta")) -or
        !(Select-String -Path ($job + ".sta") -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus H20.39 $mode solve did not report successful completion"
    }
    & "C:\SIMULIA\Commands\abaqus.bat" python extract_h20_39.py ($job + ".odb") $resultFile
    if ($LASTEXITCODE -ne 0 -or !(Test-Path $resultFile) -or (Get-Item $resultFile).Length -le 128) {
        throw "Abaqus H20.39 $mode extraction failed"
    }
    Copy-Item $resultFile $SourceDirectory
}
Write-Output "Abaqus H20.39 work directory: $work"
