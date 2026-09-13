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
$job = "b529_hex8_c3d8rt_${Geometry}_capacity_probe"
$output = "b529_hex8_c3d8rt_${Geometry}_capacity.csv"
$work = Join-Path $env:TEMP ("fuelsim_b529_" + $Geometry + "_" + (Get-Date -Format "yyyyMMdd_HHmmss"))
New-Item -ItemType Directory -Path $work | Out-Null
Copy-Item (Join-Path $SourceDirectory "${job}.inp") $work
Copy-Item (Join-Path $SourceDirectory "extract_b529.py") $work
Set-Location $work

Invoke-Abaqus job=$job input="${job}.inp" output_precision=full cpus=1 interactive
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.29 $Geometry capacity probe failed with exit code $LASTEXITCODE"
}
if (!(Test-Path "${job}.sta") -or
    !(Select-String -Path "${job}.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.29 $Geometry capacity probe did not report successful completion"
}

Invoke-Abaqus python extract_b529.py "${job}.odb" $output
if ($LASTEXITCODE -ne 0) {
    throw "Abaqus B5.29 $Geometry capacity extraction failed with exit code $LASTEXITCODE"
}
if (!(Test-Path $output)) {
    throw "Abaqus B5.29 extraction did not create $output"
}
Copy-Item $output $SourceDirectory
Write-Output "Abaqus B5.29 work directory: $work"
