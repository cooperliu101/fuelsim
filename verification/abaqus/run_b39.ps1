param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$work = Join-Path $env:TEMP ("fuelsim_b39_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_b39.py") $work
Set-Location $work

$cases = @(
    @{ Name = "matching"; NormalX = "1."; NormalY = "0." },
    @{ Name = "nonmatching"; NormalX = "1."; NormalY = "0." },
    @{ Name = "tilted_gap"; NormalX = "0.9063077870366499"; NormalY = "0.4226182617406994" }
)

foreach ($case in $cases) {
    $job = "b39_hex8_sts_" + $case.Name
    $inputFile = $job + ".inp"
    Copy-Item (Join-Path $SourceDirectory $inputFile) $work
    & "C:\SIMULIA\Commands\abaqus.bat" `
        "job=$job" `
        "input=$inputFile" `
        output_precision=full `
        interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus B3.9 $($case.Name) solve failed with exit code $LASTEXITCODE"
    }
    $statusFile = $job + ".sta"
    if (!(Test-Path $statusFile) -or
        !(Select-String -Path $statusFile -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus B3.9 $($case.Name) did not report successful completion"
    }

    $displacement = $job + "_displacement.csv"
    $force = $job + "_force.csv"
    $reaction = $job + "_reaction.csv"
    & "C:\SIMULIA\Commands\abaqus.bat" `
        python extract_b39.py `
        ($job + ".odb") `
        $inputFile `
        $case.NormalX `
        $case.NormalY `
        $displacement `
        $force `
        $reaction
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus B3.9 $($case.Name) extraction failed with exit code $LASTEXITCODE"
    }
    foreach ($result in @($displacement, $force, $reaction)) {
        if (!(Test-Path $result)) {
            throw "Abaqus B3.9 $($case.Name) extraction did not create $result"
        }
        Copy-Item $result $SourceDirectory
    }
}
Write-Output "Abaqus B3.9 work directory: $work"
