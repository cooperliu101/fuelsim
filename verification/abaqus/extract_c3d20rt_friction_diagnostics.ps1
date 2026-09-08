param([Parameter(Mandatory = $true)][string]$SourceDirectory,
      [Parameter(Mandatory = $true)][string]$WorkDirectory)
$ErrorActionPreference = 'Stop'
$Job = 'c3d20rt_finite_contact_friction'
$SourceDirectory = (Resolve-Path $SourceDirectory).Path
$WorkDirectory = (Resolve-Path $WorkDirectory).Path
if ((Get-FileHash (Join-Path $SourceDirectory "$Job.inp")).Hash -ne
    (Get-FileHash (Join-Path $WorkDirectory "$Job.inp")).Hash) { throw 'Reference input hash mismatch' }
if (!(Select-String -Path (Join-Path $WorkDirectory "$Job.sta") -Pattern 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY' -Quiet)) {
    throw 'Reference analysis did not complete'
}
Set-Location $WorkDirectory
foreach ($Item in @(@('extract_c3d20rt_contact_frames.py', 'frames'), @('extract_c3d20rt_friction_stress.py', 'shear'))) {
    Copy-Item (Join-Path $SourceDirectory $Item[0]) $WorkDirectory
    $OutputFile = "${Job}_$($Item[1]).csv"
    & 'C:\SIMULIA\Commands\abaqus.bat' python $Item[0] "$Job.odb" $OutputFile
    if ($LASTEXITCODE -ne 0) { throw 'Native contact extraction failed' }
    if ((Get-Content $OutputFile).Count -ne 3641) { throw 'Expected 3640 secondary contact samples' }
    Copy-Item $OutputFile $SourceDirectory
}
Write-Output "reference_work_directory=$WorkDirectory"
Write-Output 'matching_input=true; complete_contact_samples=3640'
