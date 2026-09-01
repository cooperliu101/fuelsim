param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$JobName = "b552_hex8_sts_cross_face"
$Work = Join-Path $env:TEMP ("fuelsim_b552_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_b522.py") $Work
Set-Location $Work
& "C:\SIMULIA\Commands\abaqus.bat" job=$JobName input="$JobName.inp" output_precision=full interactive
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.52 solve failed with exit code $LASTEXITCODE" }
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.52 did not report successful completion"
}
$Nodal = "${JobName}_nodal.csv"
$Contact = "${JobName}_contact.csv"
& "C:\SIMULIA\Commands\abaqus.bat" python extract_b522.py "$JobName.odb" $Nodal $Contact
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.52 extraction failed with exit code $LASTEXITCODE" }
foreach ($Output in @($Nodal, $Contact)) {
    if (!(Test-Path $Output)) { throw "Abaqus B5.52 extraction did not create $Output" }
    Copy-Item $Output $SourceDirectory
}
Write-Output "Completed Abaqus operator case: $JobName"
Write-Output "Abaqus B5.52 work directory: $Work"
