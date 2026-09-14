param([string]$SourceDirectory, [string]$JobName, [string]$Extractor = "extract_b13_pcmi.py")
$ErrorActionPreference = "Stop"
function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    if ($NativeExitCode -ne 0) { throw "Abaqus exited with $NativeExitCode" }
}
$Work = Join-Path $env:TEMP ("fuelsim_thermal_probe_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $SourceDirectory "*.inc") $Work
Copy-Item (Join-Path $SourceDirectory $Extractor) $Work
Write-Output "job=$JobName work_directory=$Work"
Push-Location $Work
try {
    Invoke-Abaqus job=$JobName input="$JobName.inp" output_precision=full ask_delete=OFF cpus=1 interactive
    if (!(Select-String "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) { throw "Incomplete native solve" }
    Invoke-Abaqus python $Extractor "$JobName.odb" $JobName
    foreach ($Suffix in @("_nodes.csv", "_points.csv", "_contact.csv", ".dat", ".msg", ".sta")) {
        Copy-Item "$JobName$Suffix" $SourceDirectory
    }
} finally { Pop-Location }
