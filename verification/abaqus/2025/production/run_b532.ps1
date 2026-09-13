param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory,
    [Parameter(Mandatory = $true)]
    [ValidateSet("regular", "warped")]
    [string]$Geometry
)

function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$job = "b532_hex8_c3d8rt_finite_${Geometry}_operator_probe"
$output = "b532_hex8_c3d8rt_finite_${Geometry}_operator_nodal.csv"
$integration = "b532_hex8_c3d8rt_finite_${Geometry}_operator_integration.csv"
$work = Join-Path $env:TEMP ("fuelsim_b532_" + $Geometry + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "${job}.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b532.py") $work
Set-Location $work

Invoke-Abaqus job=$job input="${job}.inp" output_precision=full cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.32 $Geometry probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${job}.sta") -or
    !(Select-String -Path "${job}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.32 $Geometry probe did not report successful completion"
}

Invoke-Abaqus python extract_b532.py "${job}.odb" $output $integration
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.32 $Geometry extraction failed with exit code $LASTEXITCODE"
}
if (!(Test-Path $output) -or !(Test-Path $integration)) {
    throw "Abaqus B5.32 $Geometry extraction did not create both reference files"
}
Copy-Item $output $SourceDirectory
Copy-Item $integration $SourceDirectory
Write-Output "Abaqus B5.32 work directory: $work"
