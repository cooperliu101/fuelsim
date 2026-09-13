$ErrorActionPreference = "Stop"
$Source = $PSScriptRoot
$Results = @()
Write-Output "BATCH_START 0 run_b7_rz.ps1"
try { & (Join-Path $Source "run_b7_rz.ps1") -Cases 'b77_rz_finite_thermal' -DestinationDirectory $Source; $Results += "0,success" } catch { Write-Output $_; $Results += "0,failed" }
Write-Output "BATCH_START 1 run_b9_rz.ps1"
try { & (Join-Path $Source "run_b9_rz.ps1") -Cases 'b91_cax4rt_finite_probe,b910_cax4rt_finite_thermal,b917_cax4rt_finite_thermal_operators' -DestinationDirectory $Source; $Results += "1,success" } catch { Write-Output $_; $Results += "1,failed" }
Write-Output "BATCH_START 2 run_b10_rz.ps1"
try { & (Join-Path $Source "run_b10_rz.ps1") -Cases 'b101_cax8t_finite_probe,b111_cax8t_finite_thermal' -DestinationDirectory $Source; $Results += "2,success" } catch { Write-Output $_; $Results += "2,failed" }
Write-Output "BATCH_START 3 run_b12_rz.ps1"
try { & (Join-Path $Source "run_b12_rz.ps1") -Cases 'b121_cax8rt_finite_probe,b1211_cax8rt_finite_thermal' -DestinationDirectory $Source; $Results += "3,success" } catch { Write-Output $_; $Results += "3,failed" }
Write-Output "BATCH_START 4 run_b15_cax4t.ps1"
try { & (Join-Path $Source "run_b15_cax4t.ps1") -Cases 'b150_cax4t_finite_thermal_operator' -DestinationDirectory $Source; $Results += "4,success" } catch { Write-Output $_; $Results += "4,failed" }
Write-Output "BATCH_START 5 run_b61.ps1"
try { & (Join-Path $Source "run_b61.ps1") -SourceDirectory $Source; $Results += "5,success" } catch { Write-Output $_; $Results += "5,failed" }
Write-Output "BATCH_START 6 run_b534.ps1"
try { & (Join-Path $Source "run_b534.ps1") -SourceDirectory $Source; $Results += "6,success" } catch { Write-Output $_; $Results += "6,failed" }
Write-Output "BATCH_START 7 run_b54.ps1"
try { & (Join-Path $Source "run_b54.ps1") -SourceDirectory $Source; $Results += "7,success" } catch { Write-Output $_; $Results += "7,failed" }
Write-Output "BATCH_START 8 run_b57.ps1"
try { & (Join-Path $Source "run_b57.ps1") -SourceDirectory $Source; $Results += "8,success" } catch { Write-Output $_; $Results += "8,failed" }
Write-Output "BATCH_START 9 run_b526.ps1"
try { & (Join-Path $Source "run_b526.ps1") -SourceDirectory $Source -Case 'b526_friction_reversal'; $Results += "9,success" } catch { Write-Output $_; $Results += "9,failed" }
Write-Output "BATCH_START 10 run_c3d20rt_probe.ps1"
try { & (Join-Path $Source "run_c3d20rt_probe.ps1") -SourceDirectory $Source -JobName 'c3d20rt_thermal_probe'; $Results += "10,success" } catch { Write-Output $_; $Results += "10,failed" }
Write-Output "BATCH_START 11 run_c3d20rt_probe.ps1"
try { & (Join-Path $Source "run_c3d20rt_probe.ps1") -SourceDirectory $Source -JobName 'c3d20rt_nonaffine_probe'; $Results += "11,success" } catch { Write-Output $_; $Results += "11,failed" }
$Results | Set-Content (Join-Path $Source "batch_status.csv")
