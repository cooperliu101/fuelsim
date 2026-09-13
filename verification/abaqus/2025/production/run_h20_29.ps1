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
$work = Join-Path $env:TEMP ("fuelsim_h20_29_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "extract_h20_29.py") $work
Set-Location $work

$cases = @(
    @{ Name = "biaxial"; NormalX = "1."; NormalY = "0." },
    @{ Name = "tilted"; NormalX = "0.9063077870366499"; NormalY = "0.4226182617406994" },
    @{ Name = "graded"; NormalX = "1."; NormalY = "0." },
    @{ Name = "soft_penalty"; NormalX = "1."; NormalY = "0." },
    @{ Name = "initial_gap"; NormalX = "1."; NormalY = "0." }
)

foreach ($case in $cases) {
    $job = "h20_29_hex20_sts_" + $case.Name
    $inputFile = $job + ".inp"
    Copy-Item (Join-Path $SourceDirectory $inputFile) $work
    Invoke-Abaqus `
        "job=$job" `
        "input=$inputFile" `
        output_precision=full `
        cpus=1 interactive
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus H20.29 $($case.Name) solve failed with exit code $LASTEXITCODE"
    }
    $statusFile = $job + ".sta"
    if (!(Test-Path $statusFile) -or
        !(Select-String -Path $statusFile -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus H20.29 $($case.Name) did not report successful completion"
    }

    $displacement = $job + "_displacement.csv"
    $force = $job + "_force.csv"
    $reaction = $job + "_reaction.csv"
    Invoke-Abaqus `
        python extract_h20_29.py `
        ($job + ".odb") `
        $inputFile `
        $case.NormalX `
        $case.NormalY `
        $displacement `
        $force `
        $reaction
    if ($LASTEXITCODE -ne 0) {
        throw "Abaqus H20.29 $($case.Name) extraction failed with exit code $LASTEXITCODE"
    }
    foreach ($result in @($displacement, $force, $reaction)) {
        if (!(Test-Path $result)) {
            throw "Abaqus H20.29 $($case.Name) extraction did not create $result"
        }
        Copy-Item $result $SourceDirectory
    }
}
Write-Output "Abaqus H20.29 work directory: $work"
