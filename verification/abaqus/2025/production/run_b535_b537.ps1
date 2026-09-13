param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDirectory,
    [string]$Case = ""
)

function Invoke-Abaqus {
    $ErrorActionPreference = "Continue"
    $NativeOutput = & "C:\SIMULIA\Commands\abq2025.bat" @args 2>&1
    $NativeExitCode = $LASTEXITCODE
    foreach ($Item in $NativeOutput) { Write-Output ($Item.ToString()) }
    $global:LASTEXITCODE = $NativeExitCode
}

$ErrorActionPreference = "Stop"
$Cases = @(
    @{ Job = "b535_hex8_c3d8rt_finite_j2"; Extractor = "extract_b510.py" },
    @{ Job = "b536_hex8_c3d8rt_finite_norton"; Extractor = "extract_b511.py" },
    @{ Job = "b537_hex8_c3d8rt_finite_coupled"; Extractor = "extract_b512.py" }
)
if ($Case -ne "") {
    $Cases = @($Cases | Where-Object { $_.Job -eq $Case })
    if ($Cases.Count -ne 1) { throw "Unknown B5.35-B5.37 case: $Case" }
}
$Work = Join-Path $env:TEMP ("fuelsim_b535_b537_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $Work | Out-Null
foreach ($Item in $Cases) {
    $JobName = $Item.Job
    $InputFile = "$JobName.inp"
    Copy-Item (Join-Path $SourceDirectory $InputFile) $Work
    Copy-Item (Join-Path $SourceDirectory $Item.Extractor) $Work
    Set-Location $Work
    Invoke-Abaqus job=$JobName input=$InputFile output_precision=full cpus=1 interactive
    if ($LASTEXITCODE -ne 0) { throw "Abaqus $JobName solve failed with exit code $LASTEXITCODE" }
    if (!(Test-Path "$JobName.sta") -or
        !(Select-String -Path "$JobName.sta" -Pattern "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" -Quiet)) {
        throw "Abaqus $JobName did not report successful completion"
    }
    $Outputs = @("${JobName}_nodal.csv", "${JobName}_integration.csv", "${JobName}_energy.csv")
    Invoke-Abaqus python $Item.Extractor "$JobName.odb" $Outputs[0] $Outputs[1] $Outputs[2]
    if ($LASTEXITCODE -ne 0) { throw "Abaqus $JobName extraction failed with exit code $LASTEXITCODE" }
    $ExpectedLines = @(81, 11, 11)
    for ($Index = 0; $Index -lt $Outputs.Count; ++$Index) {
        if (!(Test-Path $Outputs[$Index]) -or
            (Get-Content $Outputs[$Index] | Measure-Object -Line).Lines -ne $ExpectedLines[$Index]) {
            throw "Abaqus $JobName extraction produced an incomplete $($Outputs[$Index])"
        }
    }
    foreach ($Output in $Outputs) {
        Copy-Item $Output $SourceDirectory
    }
    Write-Output "Completed Abaqus material case: $JobName"
}
Write-Output "Abaqus B5.35-B5.37 work directory: $Work"
