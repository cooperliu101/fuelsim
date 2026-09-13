param(
    [string[]]$Cases = @("elastic_moduli", "elastic_moduli_finite", "plastic_yield", "plastic_yield_finite", "plastic_hardening", "plastic_hardening_finite", "creep_coefficient", "creep_coefficient_finite", "creep_exponent", "creep_exponent_finite"),
    [string]$DestinationDirectory = $PSScriptRoot
)
$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$AllowedCases = @("elastic_moduli", "elastic_moduli_finite", "plastic_yield", "plastic_yield_finite", "plastic_hardening", "plastic_hardening_finite", "creep_coefficient", "creep_coefficient_finite", "creep_exponent", "creep_exponent_finite")
New-Item -ItemType Directory -Force -Path $DestinationDirectory | Out-Null
foreach ($Case in ($Cases -join ',').Split(',')) {
    if ($Case -notin $AllowedCases) { throw "Unexpected case: $Case" }
    $Work = Join-Path $env:TEMP ("fuelsim_cax4t_material_temperature_" + [guid]::NewGuid().ToString("N"))
    New-Item -ItemType Directory -Path $Work | Out-Null
    $InputPath = Join-Path $PSScriptRoot "$Case.inp"
    Copy-Item $InputPath, (Join-Path $PSScriptRoot "extract.py") $Work
    if ((Get-FileHash $InputPath).Hash -ne (Get-FileHash (Join-Path $Work "$Case.inp")).Hash) {
        throw "Input copy changed bytes"
    }
    Write-Output "case=$Case work_directory=$Work"
    Push-Location $Work
    try {
        $ErrorActionPreference = "Continue"
        & C:\SIMULIA\Commands\abaqus.bat job=$Case input="$Case.inp" cpus=1 output_precision=full ask_delete=OFF interactive *> "${Case}_run.log"
        $AnalysisExitCode = $LASTEXITCODE
        $ErrorActionPreference = "Stop"
        foreach ($Suffix in @(".dat", ".msg", ".sta", "_run.log")) {
            if (Test-Path "$Case$Suffix") { Copy-Item "$Case$Suffix" $DestinationDirectory }
        }
        if ($AnalysisExitCode -ne 0 -or !(Test-Path "$Case.sta") -or
            !(Select-String -Path "$Case.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
            throw "Native calculation failed; inspect $Work"
        }
        $ErrorActionPreference = "Continue"
        & C:\SIMULIA\Commands\abaqus.bat python extract.py "$Case.odb" $Case *> "${Case}_extract.log"
        $ExtractExitCode = $LASTEXITCODE
        $ErrorActionPreference = "Stop"
        Copy-Item "${Case}_extract.log" $DestinationDirectory
        if ($ExtractExitCode -ne 0 -or !(Test-Path "${Case}_export_complete.txt")) {
            throw "Native extraction failed; inspect $Work"
        }
        foreach ($Suffix in @("_nodes.csv", "_points.csv", "_export_complete.txt")) {
            Copy-Item "$Case$Suffix" $DestinationDirectory
        }
        @("case=$Case", "cpus=1", "output_precision=full", "analysis_exit_code=$AnalysisExitCode",
            "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
            "input_sha256=$((Get-FileHash "$Case.inp" -Algorithm SHA256).Hash)",
            "extractor_sha256=$((Get-FileHash 'extract.py' -Algorithm SHA256).Hash)") |
            Set-Content (Join-Path $DestinationDirectory "${Case}_provenance.txt")
        Write-Output "case=$Case completed"
    } finally { Pop-Location }
}
