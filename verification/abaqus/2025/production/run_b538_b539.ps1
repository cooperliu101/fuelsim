param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory,
    [string]$Case = ""
)

function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$Work = Join-Path $env:TEMP ("fuelsim_b538_b539_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "b538_b539_cases.tsv") $Work
Copy-Item (Join-Path $SourceDirectory "extract_b524_b525.py") $Work
$Cases = Import-Csv (Join-Path $Work "b538_b539_cases.tsv") -Delimiter "`t"
if ($Case -ne "") {
    $Cases = @($Cases | Where-Object { $_.case -eq $Case })
    if ($Cases.Count -ne 1) { throw "Unknown B5.38-B5.39 case: $Case" }
}
foreach ($Item in $Cases) {
    $JobName = $Item.case
    $InputFile = "$JobName.inp"
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
    Set-Location $Work
    Invoke-Abaqus job=$JobName input=$InputFile output_precision=full cpus=1 interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus $JobName solve failed with exit code $LASTEXITCODE" }
    if (!(Test-Path "$JobName.sta") -or
        !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus $JobName did not report successful completion"
    }
    $Outputs = @("${JobName}_nodal.csv", "${JobName}_integration.csv", "${JobName}_contact.csv", "${JobName}_energy.csv")
    $Success = $Outputs[0] + ".ok"
    if (Test-Path $Success) { Remove-Item $Success }
    Invoke-Abaqus python extract_b524_b525.py "$JobName.odb" `
        $Outputs[0] $Outputs[1] $Outputs[2] $Outputs[3] $Item.expected_frames reduced
    if ($LASTEXITCODE -ne 0) { throw "Abaqus $JobName extraction failed with exit code $LASTEXITCODE" }
    if (!(Test-Path $Success)) { throw "Abaqus $JobName extraction did not reach successful completion" }
    foreach ($Output in $Outputs) {
        if (!(Test-Path $Output)) { throw "Abaqus $JobName extraction did not create $Output" }
        Copy-Item $Output $SourceDirectory
    }
    Write-Output "Completed Abaqus contact case: $JobName"
}
Write-Output "Abaqus B5.38-B5.39 work directory: $Work"
