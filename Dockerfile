FROM amazonlinux:latest AS builder

RUN yum -y install which unzip aws-cli \
	&& yum install yum-utils -y \
	&& yum-config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo \
	&& rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
	&& yum install intel-mkl-2018.2-046 -y \
	&& yum install gcc-c++ -y \
	&& yum install make -y \
	&& yum install python -y \
	&& yum install dos2unix -y \
	&& yum install python-pip -y \
	&& rm -rf /var/cache/yum/

FROM builder AS deploy	
WORKDIR /usr/src/dockerdeploy/


COPY . /usr/src/dockerdeploy/
RUN make
ENV LD_LIBRARY_PATH = $LD_LIBRARY_PATH:/opt/intel/parallel_studio_xe_2018.2.046/compilers_and_libraries_2018/linux/mkl/lib/intel64_lin/:/opt/intel/parallel_studio_xe_2018.2.046/compilers_and_libraries_2018/linux/compiler/lib/intel64_lin/:/usr/lib64/

ARG AWS_ACCESS_KEY_ID_BUILD
ARG AWS_SECRET_ACCESS_KEY_BUILD
ARG AWS_REGION_BUILD=us-east-2

ENV AWS_ACCESS_KEY_ID=${AWS_ACCESS_KEY_ID_BUILD}
ENV AWS_SECRET_ACCESS_KEY=${AWS_SECRET_ACCESS_KEY_BUILD}
ENV AWS_REGION=${AWS_REGION_BUILD}


RUN dos2unix file.sh 
RUN chmod +x file.sh
CMD ["/bin/bash", "./file.sh"]