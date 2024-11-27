module m_util
  use m_types, only: fp, i1, i2, i4
  use m_datastructures, only: SimulationParameters
  implicit none

  real(fp) :: PI = 4._fp * datan(1._fp)
  real(fp) :: k_B = 1.380649e-23_fp
  real(fp) :: R = 8.31446261815324_fp


contains

  function sigma_T(particle1, particleType1, particle2, particleType2, params) result (sigma)
    implicit none
    real(fp), dimension(5), intent(in) :: particle1, particle2
    integer(i1), intent(in) :: particleType1, particleType2
    type(SimulationParameters), intent(in) :: params

    real(fp) :: sigma, d

    d = .5_fp * (params%d_ref(particleType1)+params%d_ref(particleType2))

    if (params%collisionModel == 1) then !> hard sphere
      sigma = PI * d**2
    else
      print *, "collision_model not implemented yet: ", params%collisionModel
      sigma = -1
    end if 

  end function sigma_T


  !> calculate the intersection of two lines, where the line is given by two points: [x1,y1,x2,y2]
  !> if there is no intersection or the lines are identical, return [-1,-1]
  function lineIntersection(line1, line2) result (intersection)
    implicit none
    real(fp), dimension(4), intent(in) :: line1, line2
    real(fp), dimension(2) :: intersection
    
    real(fp) :: d, u1, u2

    !> denominator
    d = (line2(4)-line2(2))*(line1(3)-line1(1)) - (line2(3)-line2(1))*(line1(4)-line1(2))

    !> if the lines are parallel
    if (abs(d) < 1e-10_fp) then
      intersection(:) = -1
      return
    end if

    u1 = ((line2(3)-line2(1)) * (line1(2)-line2(2)) - (line2(4)-line2(2))*(line1(1)-line2(1))) / d
    u2 = ((line1(3)-line1(1)) * (line1(2)-line2(2)) - (line1(4)-line1(2))*(line1(1)-line2(1))) / d

    ! print *, "u1, u2", u1, u2
    !> if the intersection is not in the line segment defined by the points for at least one of the lines
    if (u1 < 0 .or. u1 > 1 .or. u2 < 0 .or. u2 > 1) then
      intersection(:) = -1
      return
    end if

    intersection = line1(:2) + u1 * (line1(3:4)-line1(:2))
  end function lineIntersection

  !> sort the first dimension arr1 based on arr2
  subroutine sort(arr1, arr2)
    implicit none
    real(fp), dimension(:,:), pointer, intent(inout) :: arr1
    integer(i4), dimension(:), pointer, intent(inout) :: arr2

    call quicksort(arr1, arr2)
  end subroutine sort 

  !> implementation of quick-sort where array arr1 is sorted in the first dimension based on arr2
  !> based on https://gist.github.com/1AdAstra1/6f7785373efe5bb6c254d2e20c78ccc4
  recursive subroutine quicksort(arr1, arr2)
    implicit none
    real(fp), dimension(:,:), intent(inout) :: arr1
    integer(i4), dimension(:), intent(inout) :: arr2

    real(fp), dimension(5) :: temp1
    integer(i4) :: first, last, i,j, temp2, pivot

    first = 1
    last = size(arr2,1)
    !> chose pivot element as element in the middle of the array
    pivot = arr2((first+last)/2)
    i = first
    j = last
    do
      do while(arr2(i) < pivot)
        i = i+1
      end do 
      do while(pivot < arr2(j))
        j = j-1
      end do 
      if (i >= j) then
        exit
      end if 
      temp1 = arr1(i,:)
      temp2 = arr2(i)
      arr1(i,:) = arr1(j,:)
      arr2(i) = arr2(j)
      arr1(j,:) = temp1
      arr2(j) = temp2

      i = i+1
      j = j-1
    end do

    if (first < i-1) then
      call quicksort(arr1(first:i-1,:), arr2(first:i-1))
    end if 
    if (j+1 < last) then
      call quicksort(arr1(j+1:last,:), arr2(j+1:last))
    end if 
    
  end subroutine quicksort 

end module m_util
