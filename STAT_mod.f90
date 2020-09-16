module statistics
    implicit none
    contains 
    subroutine pdf_function(N,data_sample,w_sample,pdf_points,pdf,pdf_w,pdf_b)
    ! Creates pdf using values via histogram with a input bin size over the range of values.
    implicit none
    integer(kind=4),intent(in)::N,pdf_points
    double precision,dimension(0:N-1,0:N-1,0:N-1),intent(in)::data_sample,w_sample
    double precision,dimension(0:pdf_points)::pdf,pdf_w,pdf_b
    integer(kind=4)::i_x,i_y,i_z,slot
    double precision::bin_size,data_max,data_min
    data_min=MINVAL(data_sample)
    data_max=MAXVAL(data_sample)
    bin_size=(data_max-data_min)/DBLE(pdf_points)
    pdf=0.0D0
    pdf_w=0.0D0
    pdf_b(0)=data_min
    do slot=1,pdf_points
        pdf_b(slot)=DBLE(slot-0.5D0)*bin_size+data_min
    end do
    DO i_x=0,N-1
    DO i_y=0,N-1
    DO i_z=0,N-1
        slot=FLOOR((data_sample(i_x,i_y,i_z)-data_min)/bin_size)
        pdf(slot)=pdf(slot)+1.0D0
        pdf_w(slot)=pdf_w(slot)+w_sample(i_x,i_y,i_z)
    end do
    end do
    end do
    pdf=pdf/(N**3)
    pdf_w=pdf/(N**3)
    end
    subroutine pdf_function_adv(N,data_sample,w_sample,pdf_points,data_min,data_max,pdf,pdf_w,pdf_b)
    ! Creates pdf using values via histogram with a input bin size over the range of values.
    implicit none
    integer(kind=4),intent(in)::N,pdf_points
    double precision,dimension(0:N-1,0:N-1,0:N-1),intent(in)::data_sample,w_sample
    double precision,dimension(0:pdf_points)::pdf,pdf_w,pdf_b
    double precision,intent(in)::data_max,data_min
    integer(kind=4)::i_x,i_y,i_z,slot
    double precision::bin_size
    bin_size=(data_max-data_min)/DBLE(pdf_points)
    pdf=0.0D0
    pdf_w=0.0D0
    pdf_b(0)=data_min
    do slot=1,pdf_points
        pdf_b(slot)=DBLE(slot-0.5D0)*bin_size+data_min
    end do
    DO i_x=0,N-1
    DO i_y=0,N-1
    DO i_z=0,N-1
        if (data_sample(i_x,i_y,i_z) .LT. data_min) then
            slot=0
        elseif  (data_sample(i_x,i_y,i_z) .GT. data_max) then
            slot=pdf_points
        else
            slot=FLOOR((data_sample(i_x,i_y,i_z)-data_min)/bin_size)
        end if
        pdf(slot)=pdf(slot)+1.0D0
        pdf_w(slot)=pdf_w(slot)+w_sample(i_x,i_y,i_z)**2.0
    end do
    end do
    end do
    pdf=pdf/(N**3)
    pdf_w=pdf_w/(N**3)
    end
end module statistics
