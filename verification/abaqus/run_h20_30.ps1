param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_h20_30_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_h20_30.py") $work
Set-Location $work

$cases = @("faceted_cylinder", "quadratic_cylinder")
foreach ($case in $cases) {
    $job = "h20_30_hex20_sts_" + $case
    $inputFile = $job + ".inp"
    Copy-Item (Join-Path $SourceDirectory $inputFile) $work
    & "C:\SIMULIA\Commands\abaqus.bat" `
        "job=$job" `
        "input=$inputFile" `
        output_precision=full `
        interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus H20.30 $case solve failed with exit code $LASTEXITCODE"
    }
    $statusFile = $job + ".sta"
    if (!(Test-Path $statusFile) -or
        !(Select-String -Path $statusFile -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus H20.30 $case did not report successful completion"
    }

    $displacement = $job + "_displacement.csv"
    $force = $job + "_force.csv"
    $resultant = $job + "_resultant.csv"
    & "C:\SIMULIA\Commands\abaqus.bat" `
        python extract_h20_30.py `
        ($job + ".odb") `
        $inputFile `
        $displacement `
        $force `
        $resultant
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus H20.30 $case extraction failed with exit code $LASTEXITCODE"
    }
    foreach ($result in @($displacement, $force, $resultant)) {
        if (!(Test-Path $result)) {
            throw "Abaqus H20.30 $case extraction did not create $result"
        }
        Copy-Item $result $SourceDirectory
    }
}
Write-Output "Abaqus H20.30 work directory: $work"
