module m_util
  use m_types, only: fp, i1, i2, i4
  use m_datastructures, only: SimulationParameters
  implicit none

  real(fp) :: PI = 4._fp * datan(1._fp)


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

  !> sort the second dimension arr1 based on arr2
  !> based on https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
  subroutine sort(arr1, arr2)
    implicit none
    real(fp), dimension(:,:), pointer, intent(inout) :: arr1
    integer(i4), dimension(:), pointer, intent(inout) :: arr2

    real(fp), dimension(5, (size(arr1, 2)+1)/2) :: work1
    integer(i4), dimension((size(arr1, 2)+1)/2) :: work2

    call conditionedMergeSort(arr1, arr2, work1, work2)
  end subroutine sort 

  subroutine conditionedMerge(A1, B1, C1, A2, B2, C2)
    implicit none
    real(fp), target, dimension(:, :), intent(in) :: A1, B1
    real(fp), target, dimension(:, :), intent(inout) :: C1
    integer(i4), target, dimension(:), intent(in) :: A2, B2
    integer(i4), target, dimension(:), intent(inout) :: C2

    integer(i4) :: i,j,k

    if (size(A2) + size(B2) > size(C2)) stop 1

    i = 1; j = 1
    do k = 1, size(C2)
      if (i <= size(A2) .and. j <= size(B2)) then
        if (A2(i) <= B2(j)) then
          C2(k) = A2(i)
          C1(:, k) = A1(:, i)
          i = i+1
        else
          C2(k) = B2(j)
          C1(:, k) = B1(:, j)
          j = j+1
        end if 
      else if (i <= size(A2)) then
        C1(:, k) = A1(:, i)
        C2(k) = A2(i)
        i = i+1
      else if (j <= size(B2)) then
        C1(:, k) = B1(:,j)
        C2(k) = B2(j)
        j = j+1
      end if 
    end do 
  end subroutine conditionedMerge 

  recursive subroutine conditionedMergeSort(arr1, arr2, work1, work2)
    implicit none
    real(fp), dimension(:,:), intent(inout) :: arr1, work1 
    integer(i4), dimension(:), intent(inout) :: arr2, work2

    real(fp), dimension(5) :: tmp1
    integer(i4) :: tmp2
    integer(i4) :: half
    half = (size(arr2) + 1) / 2

    if (size(arr2) < 2) then
      return
    else if (size(arr2) == 2) then
      if (arr2(1) > arr2(2)) then
        !> swap two elements
        tmp1 = arr1(:, 1)
        tmp2 = arr2(1)
        arr1(:, 1) = arr1(:, 2)
        arr2(1) = arr2(2)
        arr1(:, 2) = tmp1
        arr2(2) = tmp2
      end if
    else
      call conditionedMergeSort(arr1(:, :half), arr2(:half), work1, work2)
      call conditionedMergeSort(arr1(:, half+1:), arr2(half+1:), work1, work2)
      if (arr2(half) > arr2(half+1)) then
        work1(:, :half) = arr1(:, :half)
        work2(:half) = arr2(:half)
        call conditionedMerge(work1(:, :half), arr1(:, half+1:), arr1, work2(:half), arr2(half+1:), arr2)
      end if 
    end if 
  end subroutine conditionedMergeSort 

  !> implementation of merge-sort where array arr1 is sorted in the second dimension based on arr2
  !> both arrays are sorted inplace
  ! recursive subroutine conditionedMergeSort(arr1, arr2, left, right)
  !   real(fp), pointer, intent(inout) :: arr1(:, :)
  !   integer(i4), pointer, intent(inout) :: arr2(:)
  !   integer(i4), intent(in) :: left, right
  !   integer(i4) :: mid
  !
  !   if (left < right) then
  !     mid = (left + right) / 2
  !     call conditionedMergeSort(arr1, arr2, left, mid)
  !     call conditionedMergeSort(arr1, arr2, mid + 1, right)
  !     call conditionedMerge(arr1, arr2, left, mid, right)
  !   end if
  ! end subroutine conditionedMergeSort
  !
  ! subroutine conditionedMerge(arr1, arr2, left, mid, right)
  !   real(fp), pointer, intent(inout) :: arr1(:, :)
  !   integer(i4), pointer, intent(inout) :: arr2(:)
  !   integer(i4), intent(in) :: left, mid, right
  !   integer(i4) :: i, j, k, n1, n2, tmpMid
  !   integer(i4) :: temp1(5), temp2
  !
  !   ! n1 = mid - left + 1
  !   ! n2 = right - mid
  !
  !   tmpMid = mid
  !
  !   i = left
  !   j = mid + 1
  !   ! k = left
  !   if (arr2(mid) <= arr2(j)) then
  !     return
  !   end if 
  !
  !
  !   do while (i <= mid .and. j <= right)
  !     if (arr2(i) <= arr2(j)) then
  !       i = i + 1
  !     else
  !       temp1 = arr1(:,j)
  !       temp2 = arr2(j)
  !       do k = j, i+1, -1
  !         arr1(:, k) = arr1(:, k-1)
  !         arr2(k) = arr2(k-1)
  !       end do
  !       arr1(:, i) = temp1
  !       arr2(i) = temp2
  !       i = i + 1
  !       tmpMid = tmpMid + 1
  !       j = j + 1
  !     end if
  !   end do
  ! end subroutine conditionedMerge


end module m_util
