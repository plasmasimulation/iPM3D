module ModuleMaterialsCommunication
    use mpi
    use ModuleMaterials

    implicit none

    contains

        subroutine CommunicationMetalCharge(MT)
            type(Materials), intent(inout) :: MT
            real(8), allocatable :: charge_send(:), charge_recv(:)
            integer(4) :: i, ierr
            
            if (MT%metal_count > 0) then
                allocate(charge_send(MT%metal_count))
                allocate(charge_recv(MT%metal_count))

                do i = 1, MT%metal_count
                    charge_send(i) = MT%metls(i)%charge_one_step
                end do

                call MPI_ALLREDUCE(charge_send, charge_recv, MT%metal_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)

                do i = 1, MT%metal_count
                    MT%metls(i)%charge = MT%metls(i)%charge + charge_recv(i)
                    MT%metls(i)%charge_one_step = 0.d0
                end do

                if (allocated(charge_send)) deallocate(charge_send)
                if (allocated(charge_recv)) deallocate(charge_recv)
            end if

        end subroutine CommunicationMetalCharge


        subroutine CommunicationDielectricSigma(MT)
            type(Materials), intent(inout) :: MT
            integer(4) :: i, j, k, ierr
            real(8), allocatable :: sigma_send(:), sigma_recv(:)

            if (MT%dielectric_count > 0 .and. MT%index_length > 0) then
                allocate(sigma_send(MT%index_length))
                allocate(sigma_recv(MT%index_length))

                do k = 1, MT%index_length
                    i = MT%sigma_index(1, k)
                    j = MT%sigma_index(2, k)
                    sigma_send(k) = MT%sigma_one_step(i, j)
                end do

                call MPI_ALLREDUCE(sigma_send, sigma_recv, MT%index_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)

                do k = 1, MT%index_length
                    i = MT%sigma_index(1, k)
                    j = MT%sigma_index(2, k)
                    MT%sigma_one_step(i, j) = sigma_recv(k)
                end do

                MT%sigma = MT%sigma + MT%sigma_one_step
                MT%sigma_one_step = 0.d0

                if (allocated(sigma_send)) deallocate(sigma_send)
                if (allocated(sigma_recv)) deallocate(sigma_recv)
            end if

        end subroutine CommunicationDielectricSigma

end module ModuleMaterialsCommunication