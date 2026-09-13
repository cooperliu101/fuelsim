param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory
)

$ErrorActionPreference = "Stop"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:OPENBLAS_NUM_THREADS = "1"
$JobName = "b544_hex8_c3d8rt_distorted_bending"
$Work = Join-Path $env:TEMP ("fuelsim_b544_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory "$JobName.inp") $Work
Copy-Item (Join-Path $SourceDirectory "extract_b544.py") $Work
Set-Location $Work
& "C:\SIMULIA\Commands\abq2025.bat" job=$JobName input="$JobName.inp" cpus=1 output_precision=full interactive
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.44 solve failed with exit code $LASTEXITCODE" }
if (!(Test-Path "$JobName.sta") -or
    !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
    throw "Abaqus B5.44 did not report successful completion"
}
$Outputs = @("${JobName}_nodal.csv", "${JobName}_integration.csv", "${JobName}_contact.csv", "${JobName}_energy.csv")
$Success = $Outputs[0] + ".ok"
& "C:\SIMULIA\Commands\abq2025.bat" python extract_b544.py "$JobName.odb" `
    $Outputs[0] $Outputs[1] $Outputs[2] $Outputs[3]
if ($LASTEXITCODE -ne 0) { throw "Abaqus B5.44 extraction failed with exit code $LASTEXITCODE" }
if (!(Test-Path $Success)) { throw "Abaqus B5.44 extraction did not reach successful completion" }
foreach ($Output in $Outputs) {
    if (!(Test-Path $Output)) { throw "Abaqus B5.44 extraction did not create $Output" }
    Copy-Item $Output $SourceDirectory
}
foreach ($Extension in @("dat", "msg", "sta")) {
    Copy-Item "$JobName.$Extension" $SourceDirectory
}
@("case=$JobName", "density_input=2000 kg/m3 at initial temperature 300 K; no temperature table",
    "abaqus_version=Abaqus 2025 RELr427", "cpus=1", "output_precision=full",
    "completed_utc=$([DateTime]::UtcNow.ToString('o'))", "work_directory=$Work",
    "input_sha256=$((Get-FileHash "$JobName.inp" -Algorithm SHA256).Hash)",
    "extractor_sha256=$((Get-FileHash 'extract_b544.py' -Algorithm SHA256).Hash)") +
    @($Outputs | ForEach-Object { "${_}_sha256=$((Get-FileHash $_ -Algorithm SHA256).Hash)" }) |
    Set-Content (Join-Path $SourceDirectory "${JobName}_provenance.txt")
Write-Output "Completed Abaqus B5.44 in $Work"
