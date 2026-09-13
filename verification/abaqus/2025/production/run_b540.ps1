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
$JobName = "b540_nonmatching_contact_cycle"
$ExpectedFrames = 45
$Work = Join-Path $env:TEMP ("fuelsim_b540_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_b524_b525.py") $Work
Set-Location $Work
Invoke-Abaqus job=$JobName input="$JobName.inp" output_precision=full cpus=1 interactive
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.40 solve failed with exit code $LASTEXITCODE" }
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.40 did not report successful completion"
}
$Outputs = @("${JobName}_nodal.csv", "${JobName}_integration.csv", "${JobName}_contact.csv", "${JobName}_energy.csv")
$Success = $Outputs[0] + ".ok"
Invoke-Abaqus python extract_b524_b525.py "$JobName.odb" `
    $Outputs[0] $Outputs[1] $Outputs[2] $Outputs[3] $ExpectedFrames reduced
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.40 extraction failed with exit code $LASTEXITCODE" }
if (!(Test-Path $Success)) { throw "Abaqus B5.40 extraction did not reach successful completion" }
foreach ($Output in $Outputs) {
    if (!(Test-Path $Output)) { throw "Abaqus B5.40 extraction did not create $Output" }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Abaqus B5.40 work directory: $Work"
