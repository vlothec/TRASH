#!/usr/bin/env Rscript
print("wow")

cmd.Args = commandArgs(trailingOnly = TRUE)
change.lib.paths = TRUE
if(length(cmd.Args) > 0)
{
  if(cmd.Args[1] == "--def")
  {
    print("changed!")
  } else
  {
    print("not changed...")
  }
}