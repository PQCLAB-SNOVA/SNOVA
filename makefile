platform := avx2

subdirs :=\
  SNOVA_37_17_16_2_ESK SNOVA_37_17_16_2_SSK SNOVA_37_17_16_2_SHAKE_ESK SNOVA_37_17_16_2_SHAKE_SSK\
  SNOVA_25_8_16_3_ESK SNOVA_25_8_16_3_SSK SNOVA_25_8_16_3_SHAKE_ESK SNOVA_25_8_16_3_SHAKE_SSK\
  SNOVA_24_5_16_4_ESK SNOVA_24_5_16_4_SSK SNOVA_24_5_16_4_SHAKE_ESK SNOVA_24_5_16_4_SHAKE_SSK\
  SNOVA_56_25_16_2_ESK SNOVA_56_25_16_2_SSK SNOVA_56_25_16_2_SHAKE_ESK SNOVA_56_25_16_2_SHAKE_SSK\
  SNOVA_49_11_16_3_ESK SNOVA_49_11_16_3_SSK SNOVA_49_11_16_3_SHAKE_ESK SNOVA_49_11_16_3_SHAKE_SSK\
  SNOVA_37_8_16_4_ESK SNOVA_37_8_16_4_SSK SNOVA_37_8_16_4_SHAKE_ESK SNOVA_37_8_16_4_SHAKE_SSK\
  SNOVA_24_5_16_5_ESK SNOVA_24_5_16_5_SSK SNOVA_24_5_16_5_SHAKE_ESK SNOVA_24_5_16_5_SHAKE_SSK\
  SNOVA_75_33_16_2_ESK SNOVA_75_33_16_2_SSK SNOVA_75_33_16_2_SHAKE_ESK SNOVA_75_33_16_2_SHAKE_SSK\
  SNOVA_66_15_16_3_ESK SNOVA_66_15_16_3_SSK SNOVA_66_15_16_3_SHAKE_ESK SNOVA_66_15_16_3_SHAKE_SSK\
  SNOVA_60_10_16_4_ESK SNOVA_60_10_16_4_SSK SNOVA_60_10_16_4_SHAKE_ESK SNOVA_60_10_16_4_SHAKE_SSK\
  SNOVA_29_6_16_5_ESK SNOVA_29_6_16_5_SSK SNOVA_29_6_16_5_SHAKE_ESK SNOVA_29_6_16_5_SHAKE_SSK

OPTIMISATION = 0
valid_platform = ref

ifeq ($(platform), opt)
OPTIMISATION = 1
valid_platform = opt
endif

ifeq ($(platform), avx2)
OPTIMISATION = 2
valid_platform = avx2
endif

.PHONY: all clean $(subdirs)

all: $(subdirs)

$(subdirs): snova_config.src
	mkdir -p $@/$(valid_platform)
	grep -A 7 $@ snova_config.src > $@/$(valid_platform)/snova_config.mk
	sed -i 's/OPTIMISATION=0/OPTIMISATION=$(OPTIMISATION)/' $@/$(valid_platform)/snova_config.mk
	sh -c "cd $@/$(valid_platform) ; ln -s ../../src/* . || true"
	$(MAKE) -C $@/$(valid_platform) clean_without_libsnovasym PQCgenKAT

clean:
	rm -rf $(subdirs)
	$(MAKE) -C src/ clean_all
