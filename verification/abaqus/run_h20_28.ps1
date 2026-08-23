param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_28_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "generate_h20_28.py") $work
Copy-Item (Join-Path $SourceDirectory "extract_h20_28.py") $work
Set-Location $work

foreach ($smoothing in @("default", "linear", "quadratic")) {
    $job = "h20_28_hex20_refined_primary_sts_" + $smoothing
    $inputFile = $job + ".inp"
    $statusFile = $job + ".sta"
    $databaseFile = $job + ".odb"
    $operatorFile = $job + "_operator.csv"
    & "C:\SIMULIA\Commands\abaqus.bat" python generate_h20_28.py $inputFile $smoothing
    if (!(Test-Path $inputFile)) {
        throw "Abaqus H20.28 $smoothing input generation failed"
    }

    & "C:\SIMULIA\Commands\abaqus.bat" `
        "job=$job" `
        "input=$inputFile" `
        output_precision=full `
        interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus H20.28 $smoothing solve failed with exit code $LASTEXITCODE"
    }
    if (!(Test-Path $statusFile) -or
        !(Select-String -Path $statusFile -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus H20.28 $smoothing solve did not report successful completion"
    }

    & "C:\SIMULIA\Commands\abaqus.bat" `
        python extract_h20_28.py `
        $databaseFile `
        $operatorFile
    if (!(Test-Path $operatorFile)) {
        throw "Abaqus H20.28 $smoothing extraction did not create the operator CSV file"
    }

    Copy-Item $inputFile $SourceDirectory
    Copy-Item $operatorFile $SourceDirectory
}
Write-Output "Abaqus H20.28 work directory: $work"
