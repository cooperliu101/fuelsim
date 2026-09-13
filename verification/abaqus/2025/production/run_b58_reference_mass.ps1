param([string]$DestinationDirectory = $PSScriptRoot)
function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$JobName = "b58_hex8_c3d8t_reference_mass"
$Work = Join-Path $env:TEMP ("fuelsim_b58_reference_mass_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $PSScriptRoot "$JobName.inp"), (Join-Path $PSScriptRoot "extract_b58.py") $Work
Write-Output "case=$JobName work_directory=$Work"
Push-Location $Work
try {
    $ErrorActionPreference = "Continue"
    Invoke-Abaqus job=$JobName input="$JobName.inp" output_precision=full ask_delete=OFF cpus=1 interactive *> "${JobName}_run.log"
    $AnalysisExitCode = $LASTEXITCODE
    $ErrorActionPreference = "Stop"
    foreach ($Extension in @("dat", "msg", "sta")) {
        if (Test-Path "$JobName.$Extension") { Copy-Item "$JobName.$Extension" $DestinationDirectory }
    }
    Get-Content "${JobName}_run.log" | Set-Content (Join-Path $DestinationDirectory "${JobName}_run.log") -Encoding UTF8
    if ($AnalysisExitCode -ne 0 -or !(Test-Path "$JobName.sta") -or
        !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus failed; inspect $Work"
    }
    $ErrorActionPreference = "Continue"
    Invoke-Abaqus python extract_b58.py "$JobName.odb" "${JobName}_nodal.csv" "${JobName}_integration.csv" *> "${JobName}_extract.log"
    $ExtractExitCode = $LASTEXITCODE
    $ErrorActionPreference = "Stop"
    Get-Content "${JobName}_extract.log" | Set-Content (Join-Path $DestinationDirectory "${JobName}_extract.log") -Encoding UTF8
    if ($ExtractExitCode -ne 0 -or !(Test-Path "${JobName}_nodal.csv") -or !(Test-Path "${JobName}_integration.csv")) {
        throw "Abaqus extraction failed; inspect $Work"
    }
    Copy-Item "${JobName}_nodal.csv", "${JobName}_integration.csv" $DestinationDirectory
    @("case=$JobName", "reference_scope=initial_density_at_300K_with_unchanged_small_strain_B58_physics",
        "abaqus_version=Abaqus 2025 RELr427", "", "output_precision=full",
        "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
        "input_sha256=$((Get-FileHash "$JobName.inp" -Algorithm SHA256).Hash)",
        "extractor_sha256=$((Get-FileHash 'extract_b58.py' -Algorithm SHA256).Hash)",
        "nodal_sha256=$((Get-FileHash "${JobName}_nodal.csv" -Algorithm SHA256).Hash)",
        "integration_sha256=$((Get-FileHash "${JobName}_integration.csv" -Algorithm SHA256).Hash)") |
        Set-Content (Join-Path $DestinationDirectory "${JobName}_provenance.txt")
} finally { Pop-Location }
