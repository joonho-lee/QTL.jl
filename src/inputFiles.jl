function checkFile(file)
  #count row numbers
  f=open(file)
  nrow = countlines(f)
  close(f)

  #check column numbers
  f=open(file)
  ncol=length(split(readline(f)))
  while !eof(f)
    if length(split(readline(f))) != ncol
      error("Number of columns are not equal!!!")
    end
  end

  println("Good file with ",nrow, " rows and ",ncol," columns!!!")
  close(f)

  return nrow,ncol
end

#write out to text: writedlm is Okay, (18904,50564) matrix takes 305 seconds
#write out to binary
