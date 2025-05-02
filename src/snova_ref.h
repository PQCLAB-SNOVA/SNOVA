#ifndef SNOVA_REF_H
#define SNOVA_REF_H

/**
 * Generate private key (F part)
 * @param map2 - output: F11 F12 F21
 * @param map1 - input: P11 P12 P21 Aalpha Balpha Qalpha1 Qalpha2
 * @param T12 - input
 */
void gen_F_ref(map_group2 *map2, map_group1 *map1, T12_t T12) {
	gf16m_t temp;
	memcpy(map2->F11, map1->P11, m_SNOVA * v_SNOVA * v_SNOVA * sq_rank);
	memcpy(map2->F12, map1->P12, m_SNOVA * v_SNOVA * o_SNOVA * sq_rank);
	memcpy(map2->F21, map1->P21, m_SNOVA * o_SNOVA * v_SNOVA * sq_rank);

	for (int i = 0; i < m_SNOVA; ++i) {
		for (int j = 0; j < v_SNOVA; ++j) {
			for (int k = 0; k < o_SNOVA; ++k) {
				for (int index = 0; index < v_SNOVA; ++index) {
					gf16m_mul(map1->P11[i][j][index], T12[index][k], temp);
					gf16m_add(map2->F12[i][j][k], temp, map2->F12[i][j][k]);
				}
			}
		}
	}

	for (int i = 0; i < m_SNOVA; ++i) {
		for (int j = 0; j < o_SNOVA; ++j) {
			for (int k = 0; k < v_SNOVA; ++k) {
				for (int index = 0; index < v_SNOVA; ++index) {
					gf16m_mul(T12[index][j], map1->P11[i][index][k], temp);
					gf16m_add(map2->F21[i][j][k], temp, map2->F21[i][j][k]);
				}
			}
		}
	}

	// Clear Secret!
	SNOVA_CLEAR(temp);
}

/**
 * Generate public key (P22 part)
 * @param outP22 - output
 * @param T12 - input
 * @param P21 - input
 * @param F12 - input
 */
void gen_P22_ref(P22_byte_t outP22, T12_t T12, P21_t P21, F12_t F12) {
	gf16m_t temp1, temp2;
	P22_t P22 = {0};
	for (int i = 0; i < m_SNOVA; ++i) {
		for (int j = 0; j < o_SNOVA; ++j) {
			for (int k = 0; k < o_SNOVA; ++k) {
				for (int index = 0; index < v_SNOVA; ++index) {
					gf16m_mul(T12[index][j], F12[i][index][k], temp1);
					gf16m_mul(P21[i][j][index], T12[index][k], temp2);
					gf16m_add(temp1, temp2, temp1);
					gf16m_add(P22[i][j][k], temp1, P22[i][j][k]);
					// P22[i][j][k].printout();
				}
			}
		}
	}

	convert_GF16s_to_bytes(outP22, (uint8_t *)P22, m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA);

	// Clear Secret!
	SNOVA_CLEAR(temp1);
	SNOVA_CLEAR(temp2);
}

/**
 * evaluation_ref
 * @param hash_in_GF16Matrix - pointer to output. (Evaluation Result)
 * @param pkx - pointer to expend pk.
 * @param signature_in_GF16Matrix - pointer to signature_in_GF16Matrix.
 */
void evaluation_ref(
    gf16m_t *restrict hash_in_GF16Matrix,
    const public_key_expand *restrict pkx,
    gf16m_t *restrict signature_in_GF16Matrix
) {

	gf16m_t Left[m_SNOVA][alpha_SNOVA][n_SNOVA], Right[m_SNOVA][alpha_SNOVA][n_SNOVA];
	gf16m_t temp;

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		// evaluate signature GF16Matrix array
		for (int si = 0; si < n_SNOVA; ++si) {
			gf16m_t signature_in_GF16Matrix_transpose;
			gf16m_transpose(signature_in_GF16Matrix[si], signature_in_GF16Matrix_transpose);
			for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {

				// Left[alpha][si]
				gf16m_mul(signature_in_GF16Matrix_transpose, pkx->map1.Qalpha1[mi][alpha], temp);
				gf16m_mul(pkx->map1.Aalpha[mi][alpha], temp, Left[mi][alpha][si]);
				// Right[alpha][si]
				gf16m_mul(pkx->map1.Qalpha2[mi][alpha], signature_in_GF16Matrix[si], temp);
				gf16m_mul(temp, pkx->map1.Balpha[mi][alpha], Right[mi][alpha][si]);
			}
		}
	}

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		gf16m_set_zero(hash_in_GF16Matrix[mi]);
	}

	gf16m_t P[m_SNOVA][n_SNOVA][n_SNOVA] = {0};
	for (int mi = 0; mi < m_SNOVA; ++mi) {
		/*
		        V        O
		    +--------+--------+
		    |        |        |
		  V |  P11   |  P12   |
		    |        |        |
		    +--------+--------+   = P[n_SNOVA][n_SNOVA]
		    |        |        |
		  O |  P21   |  P22   |
		    |        |        |
		    +--------+--------+
		*/
		for (int ni = 0; ni < v_SNOVA; ++ni) {
			for (int nj = 0; nj < v_SNOVA; ++nj) {
				gf16m_clone(P[mi][ni][nj], pkx->map1.P11[mi][ni][nj]);
			}

			for (int nj = v_SNOVA; nj < n_SNOVA; ++nj) {
				gf16m_clone(P[mi][ni][nj], pkx->map1.P12[mi][ni][nj - v_SNOVA]);
			}
		}
		for (int ni = v_SNOVA; ni < n_SNOVA; ++ni) {
			for (int nj = 0; nj < v_SNOVA; ++nj) {
				gf16m_clone(P[mi][ni][nj], pkx->map1.P21[mi][ni - v_SNOVA][nj]);
			}

			for (int nj = v_SNOVA; nj < n_SNOVA; ++nj) {
				gf16m_clone(P[mi][ni][nj], pkx->P22[mi][ni - v_SNOVA][nj - v_SNOVA]);
			}
		}
	}

	for (int mi = 0; mi < m_SNOVA; ++mi) {
		// main loop
		// hash_in_GF16Matrix[i] += L[alpha][ni] * P[mi][ni]][nj] * R[alpha][nj];
		for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
			int mi_prime = i_prime(mi, alpha);

			for (int ni = 0; ni < n_SNOVA; ++ni) {
				gf16m_t sum_t0 = { 0 };
				for (int nj = 0; nj < n_SNOVA; ++nj) {
					gf16m_mul(P[mi_prime][ni][nj], Right[mi][alpha][nj], temp);
					gf16m_add(sum_t0, temp, sum_t0);
				}
				gf16m_mul(Left[mi][alpha][ni], sum_t0, temp);
				gf16m_add(hash_in_GF16Matrix[mi], temp, hash_in_GF16Matrix[mi]);
			}
		}
	}
}

/**
 * Verifies signature.
 * @param pt_digest - pointer to input digest.
 * @param pt_signature - pointer to output signature.
 * @param pk - pointer to output public key.
 * @returns - 0 if signature could be verified correctly and -1 otherwise
 */
int verify_signture_ref_core(const uint8_t *pt_digest, uint64_t bytes_digest, const uint8_t *pt_signature, const public_key_expand *pkx) {
	uint8_t hash_in_bytes[bytes_hash];
	uint8_t signed_hash[bytes_hash];
	const uint8_t *pt_salt = pt_signature + bytes_signature;

	gf16m_t hash_in_GF16Matrix[m_SNOVA];
	gf16m_t signature_in_GF16Matrix[n_SNOVA];

	Keccak_HashInstance hashInstance;
	Keccak_HashInitialize_SHAKE256(&hashInstance);
	Keccak_HashUpdate(&hashInstance, pkx->pt_public_key_seed, 8 * seed_length_public);
	Keccak_HashUpdate(&hashInstance, pt_digest, 8 * bytes_digest);
	Keccak_HashUpdate(&hashInstance, pt_salt, 8 * bytes_salt);
	Keccak_HashFinal(&hashInstance, NULL);
	Keccak_HashSqueeze(&hashInstance, signed_hash, 8 * bytes_hash);

#if (o_SNOVA * l_SNOVA) & 0x1 == 1
	signed_hash[bytes_hash - 1] &= 0x0f;
#endif

	convert_bytes_to_GF16s(pt_signature, (gf16_t *)signature_in_GF16Matrix, GF16s_signature);
	evaluation_ref(hash_in_GF16Matrix, pkx, signature_in_GF16Matrix);
	convert_GF16s_to_bytes(hash_in_bytes, (gf16_t *)hash_in_GF16Matrix, m_SNOVA * lsq_SNOVA);

	int result = 0;
	for (int i = 0; i < bytes_hash; ++i) {
		if (hash_in_bytes[i] != signed_hash[i]) {
			result = -1;
			break;
		}
	}

	return result;
}

/**
 * Regular entry point for verification
 */
int verify_signture_ref(const uint8_t *pt_digest, uint64_t bytes_digest, const uint8_t *pt_signature, const uint8_t *pk) {
	public_key_expand pkx;
	expand_public_core(&pkx, pk);
	return verify_signture_ref_core(pt_digest, bytes_digest, pt_signature, &pkx);
}

int verify_signture_pkx_ref(const uint8_t *pt_digest, uint64_t bytes_digest, const uint8_t *pt_signature, const uint8_t *pkx_pck) {
	public_key_expand pkx_unpck;
	pkx_unpack(&pkx_unpck, (public_key_expand_pack *)pkx_pck);
	return verify_signture_ref_core(pt_digest, bytes_digest, pt_signature, &pkx_unpck);
}

/**
 * Computes signature
 * @param pt_signature - pointer to output signature.
 * @param digest - pointer to input digest.
 * @param array_salt - pointer to input salt.
 * @param Aalpha -
 * @param Balpha -
 * @param Qalpha1 -
 * @param Qalpha2 -
 * @param T12 -
 * @param F11 -
 * @param F12 -
 * @param F21 -
 * @param pt_public_key_seed - pointer to output public key seed.
 * @param pt_private_key_seed - pointer to output private key seed.
 */
int sign_digest_core_ref(uint8_t *pt_signature, const uint8_t *digest, uint64_t bytes_digest, uint8_t *array_salt, Aalpha_t Aalpha,
                         Balpha_t Balpha, Qalpha1_t Qalpha1, Qalpha2_t Qalpha2, T12_t T12, F11_t F11, F12_t F12, F21_t F21,
                         const uint8_t pt_public_key_seed[seed_length_public], const uint8_t pt_private_key_seed[seed_length_private]) {
	gf16_t Gauss[m_SNOVA * lsq_SNOVA][m_SNOVA * lsq_SNOVA + 1];
	gf16_t Temp[lsq_SNOVA][lsq_SNOVA];
	gf16_t t_GF16, solution[m_SNOVA * lsq_SNOVA];

	gf16m_t Left_X_tmp, Right_X_tmp;
	gf16_t *Left_X, *Right_X;

	gf16m_t Left[m_SNOVA][alpha_SNOVA][v_SNOVA], Right[m_SNOVA][alpha_SNOVA][v_SNOVA];
	gf16m_t X_in_GF16Matrix[n_SNOVA] = {0};
	gf16m_t Fvv_in_GF16Matrix[m_SNOVA];
	gf16_t hash_in_GF16[m_SNOVA * lsq_SNOVA];
	gf16m_t signature_in_GF16Matrix[n_SNOVA];

	uint8_t signed_hash[bytes_hash];
	uint8_t vinegar_in_byte[(v_SNOVA * lsq_SNOVA + 1) >> 1];

	// temp
	gf16m_t gf16m_temp0;
	gf16m_t gf16m_temp1;
	gf16m_t gf16m_secret_temp0;

	Left_X = Left_X_tmp;
	Right_X = Right_X_tmp;
	int flag_redo = 1;
	uint8_t num_sign = 0;

	createSignedHash(digest, bytes_digest, pt_public_key_seed, array_salt, signed_hash);

	// put hash value in GF16 array
	convert_bytes_to_GF16s(signed_hash, hash_in_GF16, GF16s_hash);

	do {
		memset(Gauss, 0, sizeof(Gauss));
		num_sign++;
		flag_redo = 0;
		// put hash value in the last column of Gauss matrix
		for (int index = 0; index < (m_SNOVA * lsq_SNOVA); index++) {
			Gauss[index][m_SNOVA * lsq_SNOVA] = hash_in_GF16[index];
		}
		// generate the vinegar value
		Keccak_HashInstance hashInstance;
		Keccak_HashInitialize_SHAKE256(&hashInstance);
		Keccak_HashUpdate(&hashInstance, pt_private_key_seed, 8 * seed_length_private);
		Keccak_HashUpdate(&hashInstance, digest, 8 * bytes_digest);
		Keccak_HashUpdate(&hashInstance, array_salt, 8 * bytes_salt);
		Keccak_HashUpdate(&hashInstance, &num_sign, 8);
		Keccak_HashFinal(&hashInstance, NULL);
		Keccak_HashSqueeze(&hashInstance, vinegar_in_byte, 8 * ((v_SNOVA * lsq_SNOVA + 1) >> 1));

		convert_bytes_to_GF16s(vinegar_in_byte, (uint8_t *)X_in_GF16Matrix, v_SNOVA * lsq_SNOVA);

		// evaluate the vinegar part of central map
		for (int mi = 0; mi < m_SNOVA; ++mi) {
			for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
				for (int index = 0; index < v_SNOVA; ++index) {
					gf16m_transpose(X_in_GF16Matrix[index], gf16m_temp0);
					gf16m_mul(gf16m_temp0, Qalpha1[mi][alpha], gf16m_temp1);
					gf16m_mul(Aalpha[mi][alpha], gf16m_temp1, Left[mi][alpha][index]);
					gf16m_mul(Qalpha2[mi][alpha], X_in_GF16Matrix[index], gf16m_temp1);
					gf16m_mul(gf16m_temp1, Balpha[mi][alpha], Right[mi][alpha][index]);
				}
			}
		}
		for (int mi = 0; mi < m_SNOVA; ++mi) {
			gf16m_set_zero(Fvv_in_GF16Matrix[mi]);
		}
		for (int mi = 0; mi < m_SNOVA; ++mi) {
			for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
				int mi_prime = i_prime(mi, alpha);
				for (int j = 0; j < v_SNOVA; ++j) {
					for (int k = 0; k < v_SNOVA; ++k) {
						gf16m_mul(Left[mi][alpha][j], F11[mi_prime][j][k], gf16m_temp0);
						gf16m_mul(gf16m_temp0, Right[mi][alpha][k], gf16m_temp1);
						gf16m_add(Fvv_in_GF16Matrix[mi], gf16m_temp1, Fvv_in_GF16Matrix[mi]);
					}
				}
			}
		}
		// add to the last column of Gauss matrix
		for (int i = 0; i < m_SNOVA; ++i) {
			for (int j = 0; j < rank; ++j) {
				for (int k = 0; k < rank; ++k) {
					int index1 = i * lsq_SNOVA + j * rank + k;
					int index2 = m_SNOVA * lsq_SNOVA;
					Gauss[index1][index2] = gf16_get_add(Gauss[index1][index2], get_gf16m(Fvv_in_GF16Matrix[i], j, k));
				}
			}
		}

		// compute the coefficients of Xo and put into Gauss matrix and compute
		// the coefficients of Xo^t and add into Gauss matrix
		for (int mi = 0; mi < m_SNOVA; ++mi) {
			for (int index = 0; index < o_SNOVA; ++index) {
				for (int alpha = 0; alpha < alpha_SNOVA; ++alpha) {
					int mi_prime = i_prime(mi, alpha);
					for (int ti = 0; ti < lsq_SNOVA; ++ti) {
						for (int tj = 0; tj < lsq_SNOVA; ++tj) {
							Temp[ti][tj] = 0;
						}
					}
					for (int j = 0; j < v_SNOVA; ++j) {
						gf16m_mul(Left[mi][alpha][j], F12[mi_prime][j][index], gf16m_temp0);
						gf16m_mul(gf16m_temp0, Qalpha2[mi][alpha], Left_X_tmp);
						Left_X = Left_X_tmp;
						Right_X = Balpha[mi][alpha];
						for (int ti = 0; ti < lsq_SNOVA; ++ti) {
							for (int tj = 0; tj < lsq_SNOVA; ++tj) {
								gf16_t temp3 = 0;
								temp3 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj / rank),
								                     get_gf16m(Right_X, tj % rank, ti % rank));
								Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp3);
							}
						}
					}

					for (int j = 0; j < v_SNOVA; ++j) {
						Left_X = Aalpha[mi][alpha];
						gf16m_mul(Qalpha1[mi][alpha], F21[mi_prime][index][j], gf16m_temp0);
						gf16m_mul(gf16m_temp0, Right[mi][alpha][j], Right_X_tmp);
						Right_X = Right_X_tmp;
						for (int ti = 0; ti < lsq_SNOVA; ++ti) {
							for (int tj = 0; tj < lsq_SNOVA; ++tj) {
								gf16_t temp2 = 0;
								temp2 = gf16_get_mul(get_gf16m(Left_X, ti / rank, tj % rank),
								                     get_gf16m(Right_X, tj / rank, ti % rank));
								Temp[ti][tj] = gf16_get_add(Temp[ti][tj], temp2);
							}
						}
					}
					for (int ti = 0; ti < lsq_SNOVA; ++ti) {
						for (int tj = 0; tj < lsq_SNOVA; ++tj) {
							Gauss[mi * lsq_SNOVA + ti][index * lsq_SNOVA + tj] ^= Temp[ti][tj];
						}
					}
				}
			}
		}
		//
		// Gauss elimination
		for (int i = 0; i < m_SNOVA * lsq_SNOVA; ++i) {
			if (Gauss[i][i] == 0) {
				for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
					if (Gauss[j][i] != 0) {
						for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
							t_GF16 = Gauss[i][k];
							Gauss[i][k] = Gauss[j][k];
							Gauss[j][k] = t_GF16;
						}
						break;
					}
				}
			}
			if (Gauss[i][i] == 0) {
				flag_redo = 1;
				break;
			}

			t_GF16 = inv(Gauss[i][i]);
			for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
				Gauss[i][k] = gf16_get_mul(Gauss[i][k], t_GF16);
			}

			for (int j = i + 1; j < m_SNOVA * lsq_SNOVA; ++j) {
				if (Gauss[j][i] != 0) {
					t_GF16 = Gauss[j][i];
					for (int k = i; k < m_SNOVA * lsq_SNOVA + 1; ++k) {
						Gauss[j][k] = gf16_get_add(Gauss[j][k], gf16_get_mul(Gauss[i][k], t_GF16));
					}
				}
			}
		}

		if (!flag_redo) {
			for (int i = m_SNOVA * lsq_SNOVA - 1; i >= 0; --i) {
				t_GF16 = 0;
				for (int k = i + 1; k < m_SNOVA * lsq_SNOVA; ++k) {
					t_GF16 = gf16_get_add(t_GF16, gf16_get_mul(Gauss[i][k], solution[k]));
				}
				solution[i] = gf16_get_add(Gauss[i][m_SNOVA * lsq_SNOVA], t_GF16);
			}
		}

	} while (flag_redo);
	// printf("times of Gauss elimination : %d\n", num_sign);
	for (int index = 0; index < o_SNOVA; ++index) {
		for (int i = 0; i < rank; ++i) {
			for (int j = 0; j < rank; ++j) {
				set_gf16m(X_in_GF16Matrix[index + v_SNOVA], i, j, solution[index * lsq_SNOVA + i * rank + j]);
			}
		}
	}

	for (int index = 0; index < v_SNOVA; ++index) {
		gf16m_clone(signature_in_GF16Matrix[index], X_in_GF16Matrix[index]);
		for (int i = 0; i < o_SNOVA; ++i) {
			gf16m_mul(T12[index][i], X_in_GF16Matrix[v_SNOVA + i], gf16m_secret_temp0);
			gf16m_add(signature_in_GF16Matrix[index], gf16m_secret_temp0, signature_in_GF16Matrix[index]);
		}
	}
	for (int index = 0; index < o_SNOVA; ++index) {
		gf16m_clone(signature_in_GF16Matrix[v_SNOVA + index], X_in_GF16Matrix[v_SNOVA + index]);
	}
	// output signature
	for (int index = 0; index < n_SNOVA * lsq_SNOVA; ++index) {
		((gf16_t *)signature_in_GF16Matrix)[index] =
		    get_gf16m(signature_in_GF16Matrix[index / lsq_SNOVA], (index % lsq_SNOVA) / l_SNOVA, (index % lsq_SNOVA) % l_SNOVA);
	}
	convert_GF16s_to_bytes(pt_signature, (gf16_t *)signature_in_GF16Matrix, n_SNOVA * lsq_SNOVA);
	for (int i = 0; i < bytes_salt; ++i) {
		pt_signature[bytes_signature + i] = array_salt[i];
	}

	SNOVA_CLEAR(gf16m_secret_temp0);
	return 0;
}

#endif
