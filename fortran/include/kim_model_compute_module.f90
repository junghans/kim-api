!
! CDDL HEADER START
!
! The contents of this file are subject to the terms of the Common Development
! and Distribution License Version 1.0 (the "License").
!
! You can obtain a copy of the license at
! http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
! specific language governing permissions and limitations under the License.
!
! When distributing Covered Code, include this CDDL HEADER in each file and
! include the License file in a prominent location with the name LICENSE.CDDL.
! If applicable, add the following below this CDDL HEADER, with the fields
! enclosed by brackets "[]" replaced with your own identifying information:
!
! Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
!
! CDDL HEADER END
!

!
! Copyright (c) 2016--2018, Regents of the University of Minnesota.
! All rights reserved.
!
! Contributors:
!    Ryan S. Elliott
!

!
! Release: This file is part of the kim-api.git repository.
!


module kim_model_compute_module
  use, intrinsic :: iso_c_binding
  implicit none
  private &
    kim_model_compute_handle_equal, &
    kim_model_compute_handle_not_equal

  public &
    kim_model_compute_handle_type, &
    kim_model_compute_null_handle, &
    operator (.eq.), &
    operator (.ne.), &
    kim_model_compute_get_model_buffer_pointer, &
    kim_model_compute_log_entry, &
    kim_model_compute_string


  type, bind(c) :: kim_model_compute_handle_type
    type(c_ptr) :: p = c_null_ptr
  end type kim_model_compute_handle_type

  type, bind(c) :: kim_model_compute_type
    private
    type(c_ptr) :: p
  end type kim_model_compute_type

  type(kim_model_compute_handle_type), protected, save &
    :: kim_model_compute_null_handle

  interface operator (.eq.)
    module procedure kim_model_compute_handle_equal
  end interface operator (.eq.)

  interface operator (.ne.)
    module procedure kim_model_compute_handle_not_equal
  end interface operator (.ne.)

contains
  logical function kim_model_compute_handle_equal(left, right)
    use, intrinsic :: iso_c_binding
    implicit none
    type(kim_model_compute_handle_type), intent(in) :: left
    type(kim_model_compute_handle_type), intent(in) :: right

    if ((.not. c_associated(left%p)) .and. (.not. c_associated(left%p))) then
      kim_model_compute_handle_equal = .true.
    else
      kim_model_compute_handle_equal = c_associated(left%p, right%p)
    end if
  end function kim_model_compute_handle_equal

  logical function kim_model_compute_handle_not_equal(left, right)
    use, intrinsic :: iso_c_binding
    implicit none
    type(kim_model_compute_handle_type), intent(in) :: left
    type(kim_model_compute_handle_type), intent(in) :: right

    kim_model_compute_handle_not_equal = .not. (left .eq. right)
  end function kim_model_compute_handle_not_equal

  subroutine kim_model_compute_get_model_buffer_pointer(model_compute_handle, &
    ptr)
    use, intrinsic :: iso_c_binding
    implicit none
    interface
      subroutine get_model_buffer_pointer(model_compute, ptr) &
        bind(c, name="KIM_ModelCompute_GetModelBufferPointer")
        use, intrinsic :: iso_c_binding
        import kim_model_compute_type
        implicit none
        type(kim_model_compute_type), intent(in) :: model_compute
        type(c_ptr), intent(out) :: ptr
      end subroutine get_model_buffer_pointer
    end interface
    type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
    type(c_ptr), intent(out) :: ptr
    type(kim_model_compute_type), pointer :: model_compute

    call c_f_pointer(model_compute_handle%p, model_compute)
    call get_model_buffer_pointer(model_compute, ptr)
  end subroutine kim_model_compute_get_model_buffer_pointer

  subroutine kim_model_compute_log_entry(model_compute_handle, log_verbosity, &
    message)
    use, intrinsic :: iso_c_binding
    use kim_log_verbosity_module, only : kim_log_verbosity_type
    implicit none
    interface
      subroutine log_entry(model_compute, log_verbosity, message, line_number, &
        file_name) bind(c, name="KIM_ModelCompute_LogEntry")
        use, intrinsic :: iso_c_binding
        use kim_log_verbosity_module, only : kim_log_verbosity_type
        import kim_model_compute_type
        implicit none
        type(kim_model_compute_type), intent(in) :: model_compute
        type(kim_log_verbosity_type), intent(in), value :: log_verbosity
        character(c_char), intent(in) :: message(*)
        integer(c_int), intent(in), value :: line_number
        character(c_char), intent(in) :: file_name(*)
      end subroutine log_entry
    end interface
    type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
    type(kim_log_verbosity_type), intent(in), value :: log_verbosity
    character(len=*, kind=c_char), intent(in) :: message
    type(kim_model_compute_type), pointer :: model_compute

    call c_f_pointer(model_compute_handle%p, model_compute)
    call log_entry(model_compute, log_verbosity, trim(message)//c_null_char, &
      0, ""//c_null_char)
  end subroutine kim_model_compute_log_entry

  subroutine kim_model_compute_string(model_compute_handle, string)
    use, intrinsic :: iso_c_binding
    use kim_convert_string_module, only : kim_convert_string
    implicit none
    interface
      type(c_ptr) function model_compute_string(model_compute) &
        bind(c, name="KIM_ModelCompute_String")
        use, intrinsic :: iso_c_binding
        import kim_model_compute_type
        implicit none
        type(kim_model_compute_type), intent(in) :: model_compute
      end function model_compute_string
    end interface
    type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
    character(len=*, kind=c_char), intent(out) :: string
    type(kim_model_compute_type), pointer :: model_compute

    type(c_ptr) :: p

    call c_f_pointer(model_compute_handle%p, model_compute)
    p = model_compute_string(model_compute)
    if (c_associated(p)) then
      call kim_convert_string(p, string)
    else
      string = ""
    end if
  end subroutine kim_model_compute_string
end module kim_model_compute_module
