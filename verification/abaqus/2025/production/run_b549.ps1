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
$JobName = "b549_small_c3d20t"
$Work = Join-Path $env:TEMP ("fuelsim_b549_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($InputFile in @(
    "$JobName.inp",
    "${JobName}_mesh.inc",
    "${JobName}_mesh.json",
    "extract_b549.py"
)) {
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
}
Set-Location $Work

Invoke-Abaqus `
    job=$JobName `
    input="$JobName.inp" `
    `
    output_precision=full `
    ask_delete=OFF `
    cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.49 solve failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.49 did not report successful completion"
}

$IncrementSummary = @(
    "increment`tattempt`tsevere_discontinuity_iterations`tequilibrium_iterations`ttotal_iterations`ttime_s`tincrement_s"
)
foreach ($Line in Get-Content "$JobName.sta") {
    if ($Line -match "^\s+1\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+") {
        $Fields = ($Line.Trim() -split "\s+")
        $IncrementSummary += "$($Fields[1])`t$($Fields[2])`t$($Fields[3])`t$($Fields[4])`t$($Fields[5])`t$($Fields[7])`t$($Fields[8])"
    }
}
if ($IncrementSummary.Count -ne 21) {
    throw "Abaqus B5.49 did not produce twenty increment summary rows"
}
[IO.File]::WriteAllText(
    (Join-Path $Work "${JobName}_increments.tsv"),
    ($IncrementSummary -join "`n") + "`n",
    [Text.Encoding]::ASCII
)

Invoke-Abaqus python extract_b549.py `
    "$JobName.odb" `
    "${JobName}_mesh.json" `
    "${JobName}_temperature.csv" `
    "${JobName}_displacement.csv" `
    "${JobName}_material.csv"
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.49 extraction failed with exit code $LASTEXITCODE"
}
$ExpectedCsv = @(
    "${JobName}_temperature.csv",
    "${JobName}_displacement.csv",
    "${JobName}_material.csv"
)
foreach ($Csv in $ExpectedCsv) {
    if (!(Test-Path $Csv)) {
        throw "Abaqus B5.49 extraction did not create $Csv"
    }
}

foreach ($Output in @(
    $ExpectedCsv[0],
    $ExpectedCsv[1],
    $ExpectedCsv[2],
    "${JobName}_increments.tsv"
)) {
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B5.49 work directory: $Work"
