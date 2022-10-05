program make_mesh_for_gnuplot

character*100 :: input,output
integer :: i,j,k,i1,i2,i3
real*8,allocatable :: meshpoints(:,:)
integer,allocatable :: meshlines(:,:), meshangles(:)

call getarg(1,input)
output=trim(input)//"_triangle"

open(10,file=trim(output))

open(11,file=trim(input))
  read(11,*,err=110,end=110)
  read(11,*,err=110,end=110) i1, i2, i3
    allocate(meshpoints(i1,7)); meshpoints=0d0
    allocate(meshlines(i2,4)); meshlines=0
    allocate(meshangles(i2)); meshangles=0
  read(11,*,err=110,end=110)
  print*, i1, i2, i3
  do i=1,i1
    read(11,*,err=110,end=110) meshpoints(i,:)
  end do
  print*, "##"
  read(11,*,err=110,end=110)
  do i=1,i2
    read(11,*,err=110,end=110) meshangles(i)
    print*, i, meshangles(i)
  end do
110 close(11)
print*, "#"
open(12,file=trim(input))
  read(12,*,err=111,end=111)
  read(12,*,err=111,end=111) 
  read(12,*,err=111,end=111)
  do i=1,i1
    read(12,*,err=111,end=111) meshpoints(i,:)
  end do
  read(12,*,err=111,end=111)
  do i=1,i2
    read(12,*,err=111,end=111) k, meshlines(i,:meshangles(i))
    do j=1,k
      write(10,*) i, j, meshpoints(meshlines(i,j)+1,:)
    end do
    write(10,*) i, j, meshpoints(meshlines(i,1)+1,:)
    write(10,*) ""
  end do
111 close(12)
close(10)

end program
! splot './999out____.off_triangle' us 3:4:5 w linespoints lc "grey", './999out____.off_triangle' us 3:4:5:5 pal pt 7 w p
