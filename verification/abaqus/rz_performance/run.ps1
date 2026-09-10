param([Parameter(Mandatory=$true)][string]$SourceDirectory,
      [ValidateSet('medium','large')][string]$Size='medium',
      [switch]$Timing, [int]$Runs=1,
      [ValidateSet('cax4t','cax4rt','cax8t','cax8rt')][string]$Element='cax4t',
      [ValidateSet('small','finite')][string]$Strain='small',
      [string]$ResultsDirectory='')
$ErrorActionPreference='Stop'
if ($Element -ne 'cax4t' -and $Size -ne 'medium') { throw 'This element model has only a medium input' }
if ($Strain -eq 'finite' -and ($Size -ne 'medium' -or $Element -ne 'cax4t')) { throw 'Finite strain has only a medium CAX4T input' }
if (!$ResultsDirectory) { $ResultsDirectory = $SourceDirectory }
New-Item -ItemType Directory -Force -Path $ResultsDirectory | Out-Null
[System.Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity = [IntPtr]1
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$Work=Join-Path $env:TEMP ('fuelsim_rz_performance_'+[guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $Work | Out-Null
Copy-Item (Join-Path $SourceDirectory '*.inc') $Work
Copy-Item (Join-Path $SourceDirectory '*.inp') $Work
Copy-Item (Join-Path $SourceDirectory 'extract_results.py') $Work
$Job='rz_performance_'+$Size
if ($Element -ne 'cax4t') { $Job += '_'+$Element }
if ($Strain -eq 'finite') { $Job += '_cax4t_finite' }
if ($Timing) { $Job += '_timing' }
Write-Output "work_directory=$Work"
Push-Location $Work
try {
  for ($Run=0; $Run -lt $Runs; $Run++) {
    foreach ($Ext in @('sta','dat','msg','odb','sim','prt','com','log')) {
      $Previous = Join-Path $Work "$Job.$Ext"
      if (Test-Path $Previous) { Remove-Item $Previous }
    }
    $Clock=[System.Diagnostics.Stopwatch]::StartNew()
    & C:\SIMULIA\Commands\abaqus.bat job=$Job input="$Job.inp" cpus=1 output_precision=full ask_delete=OFF interactive
    $Code=$LASTEXITCODE
    $Clock.Stop()
    if ($Code -ne 0 -or !(Test-Path "$Job.sta") -or !(Select-String "$Job.sta" -Pattern 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' -Quiet)) { throw "Abaqus failed in $Work" }
    @("run=$Run", "external_wall_seconds=$($Clock.Elapsed.TotalSeconds.ToString('F6',[Globalization.CultureInfo]::InvariantCulture))", "cpus=1", "affinity_mask=1", "work_directory=$Work") | Set-Content (Join-Path $ResultsDirectory "${Job}_run${Run}_timing.txt")
    foreach ($Ext in @('dat','msg','sta')) { Copy-Item "$Job.$Ext" (Join-Path $ResultsDirectory "${Job}_run${Run}.$Ext") }
  }
  if (!$Timing) {
    & C:\SIMULIA\Commands\abaqus.bat python extract_results.py "$Job.odb" $Job
    if ($LASTEXITCODE -ne 0) { throw 'Extraction failed' }
    foreach ($Kind in @('nodes','points','contact')) { Copy-Item "${Job}_${Kind}.csv.gz" $ResultsDirectory }
  }
} finally { Pop-Location }
