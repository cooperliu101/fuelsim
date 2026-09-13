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
$JobName = "b552_hex8_sts_cross_face"
$Work = Join-Path $env:TEMP ("fuelsim_b552_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_b522.py") $Work
Set-Location $Work
Invoke-Abaqus job=$JobName input="$JobName.inp" output_precision=full cpus=1 interactive
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.52 solve failed with exit code $LASTEXITCODE" }
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.52 did not report successful completion"
}
$Nodal = "${JobName}_nodal.csv"
$Contact = "${JobName}_contact.csv"
Invoke-Abaqus python extract_b522.py "$JobName.odb" $Nodal $Contact
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.52 extraction failed with exit code $LASTEXITCODE" }
foreach ($Output in @($Nodal, $Contact)) {
    if (!(Test-Path $Output)) { throw "Abaqus B5.52 extraction did not create $Output" }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Completed Abaqus operator case: $JobName"
Write-Output "Abaqus B5.52 work directory: $Work"
