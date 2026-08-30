param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b544_hex8_c3d8rt_distorted_bending"
$Work = Join-Path $env:TEMP ("fuelsim_b544_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_b544.py") $Work
Set-Location $Work
& "C:\SIMULIA\Commands\abaqus.bat" job=$JobName input="$JobName.inp" output_precision=full interactive
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.44 solve failed with exit code $LASTEXITCODE" }
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.44 did not report successful completion"
}
$Outputs = @("${JobName}_nodal.csv", "${JobName}_integration.csv", "${JobName}_contact.csv", "${JobName}_energy.csv")
$Success = $Outputs[0] + ".ok"
& "C:\SIMULIA\Commands\abaqus.bat" python extract_b544.py "$JobName.odb" `
    $Outputs[0] $Outputs[1] $Outputs[2] $Outputs[3]
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.44 extraction failed with exit code $LASTEXITCODE" }
if (!(Test-Path $Success)) { throw "Abaqus B5.44 extraction did not reach successful completion" }
foreach ($Output in $Outputs) {
    if (!(Test-Path $Output)) { throw "Abaqus B5.44 extraction did not create $Output" }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Completed Abaqus B5.44 in $Work"
